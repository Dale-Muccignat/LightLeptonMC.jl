module Structs

using AutoParameters, Gases, StaticArrays, AndExport, Interpolations, QuadGK, MCSimUnits

export MCMODEL, MCPARAMS, PART, MEAS, MEAST, MEASS, MEASC, findMinInd

"""
NonConservativeMethod

Parent type of each option for ways to handle non-conservative collisions

"""
@xport abstract type NonConservativeMethod end
@xport struct Weighting <: NonConservativeMethod end
@xport struct Explicit <: NonConservativeMethod end

"""
Velocity

Parent type for each option for velocity. See subtypes for a list of options
```julia
    ZeroVelocity()
    InitVelocity()
    GaussianEnergyVelocity()
    GaussianDirectionVelocity()
```

"""
@xport abstract type Velocity end

"""
Spatial

Parent type for each option for spatial config. See subtypes for a list of options
```julia
    ZeroSpatial()
    InitSpatial(val,axis)
    GaussianLμJSpatial(val, Δy_source,gaussian_axis,axis)
```

"""
@xport abstract type Spatial end 
@xport abstract type Temporal end 
@xport abstract type Density end
@xport abstract type Field end
@xport abstract type Temp end
@xport abstract type Ionisation end

@xport struct ZeroSpatial <: Spatial end

@xport @AutoParm struct InitSpatial <: Spatial
    val::AUTO
    axis::AUTO
end

@xport @AutoParm struct InitTemporal <: Temporal
    val::AUTO = 0. *s
end

@xport @AutoParm struct GaussianTemporal <: Spatial
    val::AUTO = 20. *ns
    Δt_source::AUTO = 10. *ns
end

@xport @AutoParm struct GaussianLμJSpatial <: Spatial # TODO: Rename to remove LuJ
    val::AUTO = -7.6e-6m
    Δy_source::AUTO = 7.5e-6m
    gaussian_axis::AUTO = [0,1,0]
    axis::AUTO = [0,0,1]
end

@xport function GaussianLμJSpatial(;val=-7.5μm,Δy_source=7.5μm,gaussian_axis=[0,1,0],axis=[0,0,1])
    gaussian_axis = @SVector [gaussian_axis[1],gaussian_axis[2],gaussian_axis[3]]
    axis = @SVector [axis[1],axis[2],axis[3]]
    return GaussianLμJSpatial(val,Δy_source,gaussian_axis,axis)
end

@xport @AutoParm struct FunctionalSpatial <: Spatial
    getAxis::AUTO
    getPosition::AUTO
end

"""
FunctionalSpatial(spatialDist;axis)

spatialDist - Two element tuple. First element is the dependent variable of the density vector (second element)
axis - function that returns the axis for the distribution, can be a random vector to simulate symmetry
"""
@xport function FunctionalSpatial(spatialDist::Tuple{Any,Any},axis=()->[0,0,1])
    # Takes a two element vector spatial histogram/density distribution, creates a cumulative distribution function
    # assumes isotropic scattering
    dependent = spatialDist[1]
    independent = spatialDist[2]
    independent_func = LinearInterpolation(dependent,independent)
    vals = zero(dependent)

    norm = quadgk(independent_func,dependent[1],dependent[end],order=100)[1]

    # vals[1] = quadgk(independent,dependent[1],dependent[1],order=100)[1]
    vals[1] = zero(independent[1]*dependent[1])
    for i in 2:(length(dependent)-1)
        vals[i] = quadgk(independent_func,dependent[i-1],dependent[i],order=100)[1]
    end
    vals = vals ./ norm
    vals[end] = 1. - sum(vals)
    vals = cumsum(vals)

    spl(x) = LinearInterpolation(Interpolations.deduplicate_knots!(vals),dependent,extrapolation_bc=Flat())(x)

    return FunctionalSpatial(axis,spl)
end

@xport @AutoParm struct UniformDensity <: Density
    val::AUTO = 1e23m^-3
end
@xport @AutoParm struct LμJDensityGradient <: Density
    val::AUTO = 1e28m^-3
    background::AUTO = 2e18m^-3
    func::AUTO = (x) -> 0.5*(1e28m^-3+1e23m^-3)-0.5*(1e28m^-3-1e23m^-3)*tanh(2*ustrip(x)/1e-9)
    δtol::AUTO = 0.1
    δphase::AUTO = 1nm
end
@xport function LμJDensityGradient(;val=1e28m^-3,background=2e18m^-3,jet_width=15μm,δphase=2nm,δtol=0.01)
    func = (x) -> 0.5*(val+background)-0.5*(val-background)*tanh(2*(x-jet_width)/2/δphase)
    return LμJDensityGradient(val,background,func,δtol,δphase)
end
@xport @AutoParm struct LμJDensity <: Density
    val::AUTO = 1e28m^-3
    background::AUTO = 2e18m^-3
end
@xport @AutoParm struct CellDensity <: Density
    val::AUTO = 1e28m^-3
end
@xport @AutoParm struct ZeroField <: Field
    val::AUTO = 0. *Td
end
@xport @AutoParm struct UniformField <: Field
    val::AUTO = 1. *Td
end
@xport @AutoParm struct UniformTemp <: Temp
    val::AUTO = 0. *K
end

@xport @AutoParm struct UniformIonisation <: Ionisation 
    method::Symbol = :WEIGHTING
end

@xport @AutoParm struct RatioIonisation <: Ionisation
    method::Symbol = :WEIGHTING
    ratio::AUTO = (e) -> 0.5
end

@xport @AutoParm struct FunctionalIonisation <: Ionisation
    method::Symbol = :WEIGHTING
    func::AUTO = (e) -> 0.5
end

@xport function FunctionalIonisationInit(;func=(e,E2) -> 1.,thresh=10.,method=:WEIGHTING)
    cbins = [0.,(10. .^(-3:0.1:6.))...] .*eV #TODO: This is a dubplicate and depends on the initParams funciton
    cbins2 = [0.,(10. .^(-3:0.1:6.))...] .*eV
    matrix::Matrix{Float64} = zeros(length(cbins),length(cbins2))
    for i in eachindex(cbins)
        func2(E2) = func(cbins[i],E2) < 1e-10 ? 1e-10 : func(cbins[i],E2)
        if cbins[i] > thresh
            norm = quadgk(func2,0,cbins[i]-thresh)[1]
        else
            norm = 1
        end

        for j in eachindex(cbins2)
            if j == 1
                matrix[i,j] = quadgk(func2,0,cbins2[j])[1]/norm
            elseif j == length(cbins2)
                # matrix[i,j] = 1. - sum(matrix[i,:])
            else
                matrix[i,j] = quadgk(func2,cbins2[j-1],cbins2[j])[1]/norm
            end
        end
    end
    # matrix = replace_0(replace_nan(matrix))
    matrix2 = cumsum(matrix,dims=2)
    spl2 = LinearInterpolation((cbins,cbins2),matrix2,extrapolation_bc=Flat())	#2D spline of the matrix
    spl3 = [LinearInterpolation(Interpolations.deduplicate_knots!(spl2(c,cbins2),move_knots=true),cbins2,extrapolation_bc=Flat()) for c in cbins]
    return FunctionalIonisation(spl3,method)#, matrix,matrix2
end


@xport @AutoParm struct DiscreteIonisation <: Ionisation
    method::Symbol = :WEIGHTING
    loss::AUTO = 5. *eV
end

@xport @AutoParm struct DiscreteFunctionalIonisation <: Ionisation
    method::Symbol = :WEIGHTING
    loss::AUTO = (e) -> 5. *eV
end

@xport struct ZeroVelocity <: Velocity end

@xport @AutoParm struct InitVelocity <: Velocity
    ε::AUTO = 100. *eV
    axis::AUTO = [0,0,1]
end

@xport @AutoParm struct GaussianEnergyVelocity <: Velocity
    ε::AUTO = 10. *eV
    Δε_source::AUTO = 0.5eV
    axis::AUTO = [0,0,1]
end

@xport @AutoParm struct GaussianDirectionVelocity <: Velocity
    ε::AUTO = 10. *eV
end

@xport @AutoParm struct FunctionalEnergy <: Velocity
    getAxis::AUTO
    getVelocity::AUTO
end

@xport function FunctionalEnergy(velocityDist)
    # takes velocityDist and creates a cumulative distribution function to sample from
    # assumes isotropic scattering
    energies = 10. .^(-4:0.001:2) .*sqrt(eV)
    # energies = velocityDist.itp.knots[1]
    # P = velocityDist.itp.coefs
    # energies = sqrt.(0eV:1e-5eV:100eV)
    vals = zeros(length(energies)) 

    norm = quadgk(velocityDist,sqrt(0eV),energies[end],order=100)[1]
    # norm = energies .* P |> sum
    # @show norm
    # vals = energies .* P
    # vals ./= sum(vals)

    vals[1] = quadgk(velocityDist,sqrt(0eV),energies[1],order=100)[1]
    for i in 2:(length(energies)-1)
        vals[i] = quadgk(velocityDist,energies[i-1],energies[i],order=100)[1]
        vals[i] < 0 && (vals[i] = 0)
    end
    # return vals
    @show sum(vals)
    vals = vals ./ norm
    @show sum(vals)
    vals = cumsum(vals)
    vals[end] = 1. 

    function iso() # TODO: Moves these outside of the function
        coschi = 1-2*rand()
        sinchi = sqrt(1-coschi^2)
        phi = rand()*2*pi
        return [cos(phi)*sinchi,sin(phi)*sinchi,coschi]
    end
    function delta()
        return [0,0,1]
    end
    @show energies[1], vals[1]

    spl = LinearInterpolation(Interpolations.deduplicate_knots!(vals),energies.^2,extrapolation_bc=0eV)

    # return FunctionalEnergy(delta,spl)
    return FunctionalEnergy(iso,spl)
end

"""
MCMODEL

MCMODEL is the main constructor for MCSim.

```julia
    type = SWARM()                   # Used to set up specific models: LμJ(), SWARM() or CELL().
    particle_type = :electron        # Change particle type: :electron or :positron. This will modify ionisation behaviour.
    n = 100                          # Number of particles to simulate
    v_0 = ZeroVelocity()             # Initial velocity: ZeroVelocity()...
    r_0 = ZeroSpatial()              # Initial spatial: ZeroSpatial()...
    n0 = UniformDensity()            # Density distribution:
    En0 = UniformField()             # Electric field distribution:
    T = UniformTemp()                # Temperature distribution:
    ion_type = UniformIonisation()   # Ionisation sharing type:
    meas = [:t]                      # Variables to measure. There are many possibilities, see ?MEASLμJ, ?MEASCELL or ?MEAST
```

"""
@AutoParm mutable struct MCMODEL
    type::AUTO = SWARM()
    particle_type::AUTO = :electron
    n::Int64 = 100
    v_0::AUTO = ZeroVelocity()
    r_0::AUTO = ZeroSpatial()
    t_0::AUTO = InitTemporal()
    n0::AUTO = UniformDensity()
    En0::AUTO = UniformField()
    T::AUTO = UniformTemp()
    ion_type::AUTO = UniformIonisation()
    non_conservative_method::AUTO = Weighting()
    meas::AUTO = [:t]
    optimise_spatial::AUTO = true
end

@xport @AutoParm mutable struct SWARM
    t_f::AUTO = 2e-5s
end

"""
LμJ

This struct is used to initilize a micro-jet simulation. Used in `MCMODEL(type=LμJ())`.

Parameters:
```julia
    Δε_detector = 0.9       # Energy bin width for the detector
    Δχ_detector = 2         # Angle bin width for the detector in degrees
    thresh = 1.             # Set the minimum energy to simulate electrons too
    jet_width = 15e-6       # Diameter for the liquid jet
    jet_center = [0.,0.,0.] # Set's an offset for the jet_center
```

For more information, type ?MCMODEL
"""
@xport @AutoParm mutable struct LμJ
    Δε_detector::AUTO = 0.9eV
    Δχ_detector::AUTO = 2°
    thresh::AUTO = 1. *eV
    jet_width::AUTO = 15e-6m
    jet_center::AUTO = [0.,0.,0.] .*m
end

"""
CELL

Defines the parameters for a cell (rectangular prism) style experiment
```julia
    Δε_detector = 0.9  # energy resolution of the detector
    thresh = 1.        # minimum energy to simulate
    z = 15e-6          # depth of the cell
    z_n = 100          # number of z bins between 0 and z
    x = 15e-6          # height of cell
    y = 15e-6          # width of cell
    bin_no = 100       # no_bins for radial measurement along back surface
    r = 1e-6           # radius for measurements along the back surface of the cell
```
"""
@xport @AutoParm mutable struct CELL
    # Defined spatial parameters for the CELL
    Δε_detector::AUTO = 0.9eV
    thresh::AUTO = 1. *eV
    z::AUTO = 15e-6m
    z_n::AUTO = 100
    x::AUTO = 15e-6m
    y::AUTO = 15e-6m
    bin_no::AUTO = 100
    r::AUTO = 1e-6m
    t_f::AUTO = 2e-5s
end

@xport function LμJ(;jet_center=[0.,0.,0.].*m,Δε_detector=0.9eV,Δχ_detector=2°,thresh=1. *eV,jet_width=15μm)
    LμJ(Δε_detector,Δχ_detector,thresh,jet_width,jet_center)
end

@xport function CELL(;Δε_detector=0.9eV,thresh=1. *eV,z=15μm,z_n=100,x=15μm,y=15μm,bin_no=100,r=1μm,t_f=2e-5s)
    CELL(Δε_detector,thresh,z,z_n,x,y,bin_no,r,t_f)
end

@AutoParm mutable struct MCPARAMS
    cbins::AUTO = [0.,(10. .^(-3:0.01:2))...] .*eV
    ebins::AUTO = nothing
    gas::AUTO = nothing
    splat_gas::AUTO = nothing
    total_freq::AUTO = nothing
    a::AUTO = nothing
    # n0::Union{Vector{Float64},Float64,Nothing} = nothing
    n0::AUTO = nothing
    rbins::AUTO = nothing
    χbins::AUTO = nothing
    tbins::AUTO = collect(LinRange(0,1e-5,10)) .*s
    zbins::AUTO = nothing
    xbins::AUTO = nothing
    ybins::AUTO = nothing
    length_scale::AUTO = nothing
    time_scale::AUTO = nothing
    energy_scale::AUTO = nothing

end

function MCPARAMS(;cbins=nothing,ebins=nothing,gas=nothing,splat_gas=nothing,total_freq=nothing,a=nothing,n0=nothing,rbins=nothing,χbins=nothing,tbins=nothing,zbins=nothing,xbins=nothing,ybins=nothing,length_scale=nothing,time_scale=nothing,energy_scale=nothing)
    MCPARAMS(cbins,ebins,gas,splat_gas,total_freq,a,n0,rbins,χbins,tbins,zbins,xbins,ybins,length_scale,time_scale,energy_scale)
end


@xport @AutoParm mutable struct PART
    type::Symbol = :electron
    alive::Bool = true
    v::AUTO = [0.,0.,0.]
    r::AUTO = [0.,0.,0.]
    t::AUTO = 0. *s̃
    weight::AUTO = 1.
    tbin::Int64 = 1
    cbin::Int64 = 1
    ebin::Int64 = 1
    zbin::Int64 = 1
    ybin::Int64 = 1
    xbin::Int64 = 1
    rbin::Int64 = 1
end

@xport function PART()

    r = SArray{Tuple{3},Float64,1,3}(0.,0.,0.)
    r = r*m̃
    v = SArray{Tuple{3},Float64,1,3}(0.,0.,0.)
    v = v*m̃/s̃
    return PART(r=r,v=v)
end


@AutoParm mutable struct MEAS
    measLμJ::AUTO = MEASLμJ()
    measCell::AUTO = MEASCELL()
    measg::AUTO = MEASG()
    meast::AUTO = MEAST()
    meass::AUTO = MEASS()
    measc::AUTO = MEASC()
end

function MEAS(type::SWARM,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(),MEASCELL(),MEASG(model,params),MEAST(params),MEASS(),MEASC(params))
end

function MEAS(type::LμJ,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(model,params),MEASCELL(),MEASG(),MEAST(),MEASS(),MEASC())
end

function MEAS(type::CELL,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(),MEASCELL(model,params),MEASG(model,params),MEAST(params),MEASS(),MEASC(params))
end

@xport @AutoParm mutable struct MEASLμJ
    ε::AUTO = nothing
    χ::AUTO = nothing
    εχ::AUTO = nothing
    εχ_reject::AUTO = nothing
    zyt::AUTO = nothing
end

function MEASLμJ(model::MCMODEL,params::MCPARAMS)
    ε,χ,zyt,εχ,εχ_reject = nothing, nothing, nothing, nothing, nothing
    if :zyt ∈ model.meas
        zyt= zeros(length(params.zbins),length(params.ybins),length(params.tbins))
    end

    if :ε ∈ model.meas
        ε = zeros(length(params.ebins))
    end

    if :χ ∈ model.meas
        χ = zeros(length(params.χbins))
    end

    if :εχ ∈ model.meas || :εχ_ ∈ model.meas
        εχ = zeros(length(params.ebins),length(params.χbins))
        εχ_reject = zeros(length(params.ebins),length(params.χbins))
    end

    return MEASLμJ(ε, χ, εχ,εχ_reject,zyt)

end

@xport @AutoParm mutable struct MEASCELL
    ε::AUTO = nothing
    εz::AUTO = nothing
    εz_w::AUTO = nothing
    x::AUTO = nothing
    y::AUTO = nothing
    z::AUTO = nothing
    xy::AUTO = nothing
    xybins::AUTO = nothing
    tz::AUTO = nothing
    v_tz::AUTO = nothing
end

function MEASCELL(model::MCMODEL,params::MCPARAMS)
    εz_w, ε, εz, x, y, z, xy, xybins, tz, v_tz = nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    if :ε ∈ model.meas
        ε = zeros(length(params.ebins))
    end
    if :εz ∈ model.meas
        εz = zeros(length(params.zbins)).*eV
        εz_w = zeros(length(params.zbins))
    end
    if :x ∈ model.meas
        x = zeros(length(params.xbins))
    end
    if :y ∈ model.meas
        y = zeros(length(params.ybins))
    end
    if :z ∈ model.meas
        z = zeros(length(params.zbins))
    end
    if :xy ∈ model.meas
        xy = zeros(model.type.bin_no)
        xybins = convert(Vector{Float64},LinRange(0,model.type.r,model.type.bin_no))
    end
    if :tz ∈ model.meas
        tz = zeros(length(params.tbins),length(params.zbins))
        v_tz = zeros(length(params.tbins),length(params.zbins))*m/s
    end

    return MEASCELL(ε,εz,εz_w,x,y,z,xy,xybins,tz,v_tz)
end

@AutoParm mutable struct MEASG
    coschi::AUTO = 0.
    coll_count::AUTO = []
    SP_ε::AUTO = nothing
    SP_r::AUTO = nothing
    SP_εs::AUTO = nothing
end

function MEASG(model::MCMODEL,params::MCPARAMS)
    coschi, coll_count, SP_ε, SP_r, SP_εs, SP_zs = 0., zeros(length(params.splat_gas.cs)), nothing, nothing, nothing, nothing
    @debug coll_count
    if :SP ∈ model.meas
        SP_ε = 0.
        SP_r = 0.
    end
    if :SP_z ∈ model.meas
        SP_εs = zeros(length(params.zbins))
    end
    return MEASG(coschi, coll_count, SP_ε, SP_r, SP_εs)
end

@AutoParm mutable struct MEAST
    εᵢ::AUTO = nothing
    r::AUTO = nothing
    v::AUTO = nothing
    r²::AUTO = nothing
    v²::AUTO = nothing
    rv::AUTO = nothing
    w::AUTO = nothing
    zeff::AUTO = nothing
end

function MEAST(params)

    MEAST(zeros(length(params.tbins)-1).*eV,
          zeros(length(params.tbins)-1,3).*m,
          zeros(length(params.tbins)-1,3).*m/s,
          zeros(length(params.tbins)-1,3).*m^2,
          zeros(length(params.tbins)-1,3).*m^2/s^2,
          zeros(length(params.tbins)-1,3).*m^2/s,
          zeros(length(params.tbins)-1),
          zeros(length(params.tbins)-1)
          )

end

@AutoParm mutable struct MEASS

end

@AutoParm mutable struct MEASC
    εᵢ_final::AUTO = nothing
    r_final::AUTO = nothing
    r²_final::AUTO = nothing
    v_final::AUTO = nothing
    rv_final::AUTO = nothing
    r²v_final::AUTO = nothing
    rx²vz_final::AUTO = nothing
    rzrxvx_final::AUTO = nothing
    εᵢ::AUTO = nothing
    r::AUTO = nothing
    r²::AUTO = nothing
    v::AUTO = nothing
    rv::AUTO = nothing
    r²v::AUTO = nothing
    rx²vz::AUTO = nothing
    rzrxvx::AUTO = nothing
    weight::AUTO = nothing
    dt::AUTO = nothing
end

function MEASC(params)
    return MEASC(
          zeros(length(params.tbins)).*eV,
          zeros(length(params.tbins),3).*m,
          zeros(length(params.tbins),3).*m.^2,
          zeros(length(params.tbins),3).*m./s,
          zeros(length(params.tbins),3).*m^2 ./s,
          zeros(length(params.tbins),3).*m^3 ./s,
          zeros(length(params.tbins)).*m^3 ./s,
          zeros(length(params.tbins)).*m^3 ./s,
          zeros(length(params.tbins)).*eV*s,
          zeros(length(params.tbins),3).*m*s,
          zeros(length(params.tbins),3).*m^2*s,
          zeros(length(params.tbins),3).*m,
          zeros(length(params.tbins),3).*m^2,
          zeros(length(params.tbins),3).*m^3,
          zeros(length(params.tbins)).*m^3,
          zeros(length(params.tbins)).*m^3,
          zeros(length(params.tbins)).*s,
          zeros(length(params.tbins)).*s
          )
end


@AutoParm @xport mutable struct SWARMMEAS
    ε::AUTO = 0.
    εs::AUTO = [0.]
    R::AUTO = 0.
    Rs::AUTO = [0.]
    bW::AUTO = 0.
    bWs::AUTO = [0.]
    fW::AUTO = 0.
    fWs::AUTO = [0.]
    bDt::AUTO = 0.
    bDts::AUTO = [0.]
    fDt::AUTO = 0.
    fDts::AUTO = [0.]
    bDl::AUTO = 0.
    bDls::AUTO = [0.]
    fDl::AUTO = 0.
    fDls::AUTO = [0.]
    fQl::AUTO = 0.
    fQls::AUTO = [0.]
    fQt::AUTO = 0.
    fQts::AUTO = [0.]
end

@fastmath function findMinInd(x::AbstractArray,y)::Int64
    left::Int64 = 1
    right::Int64 = length(x)
    x = collect(x)
    d= Inf*x[1]
    dlast= -Inf*x[1]
	c::Int64 = 0
	while d != dlast
		# c = Int(round((left+right)/2))
		c = div(left+right,2)
		dlast = d
		@inbounds d = x[c]
		if d < y
			left = c
        elseif d >= y
            right = c
        end

    end
    if x[right] <= y
        return right
    else
        return left
    end
end

@fastmath function findMinInd(x::Array{Int64,1},y::Int64)::Int64
    left::Int64 = 1
    right::Int64 = length(x)
    d::Float64  = Inf
    dlast::Float64  = -Inf
    c::Int64 = 0
    while d != dlast
        # c = Int(round((left+right)/2))
        c = div(left+right,2)
        dlast = d
        @inbounds d = x[c]
        if d < y
            left = c
        elseif d >= y
            right = c
        end

    end
    if x[right] <= y
        return right
    else
        return left
    end
end
end # module
