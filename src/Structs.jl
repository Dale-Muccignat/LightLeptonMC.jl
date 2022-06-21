module Structs

using AutoParameters, Gases, StaticArrays, AndExport, Interpolations, QuadGK, Unitful, Unitful.DefaultSymbols

export MCMODEL, MCINIT, MCPARAMS, PART, MEAS, MEAST, MEASS, MEASC, findMinInd

"""
Velocity

Parent type for each option for velocity. See subtypes for a list of options
"""
@xport abstract type Velocity end
@xport abstract type Spacial end
@xport abstract type Density end
@xport abstract type Field end
@xport abstract type Temp end
@xport abstract type Ionisation end

@xport struct ZeroSpacial <: Spacial end
@xport @AutoParm struct InitSpacial <: Spacial
    val::AUTO
    axis::AUTO
end
@xport @AutoParm struct GaussianLμJSpacial <: Spacial
    val::AUTO = -7.6e-6m
    Δy_source::AUTO = 7.5e-6m
    gaussian_axis::AUTO = [0,1,0]
    axis::AUTO = [0,0,1]
end

@xport @AutoParm struct UniformDensity <: Density
    val::AUTO = 1e23
end
@xport @AutoParm struct LμJDensity <: Density
    val::AUTO = 1e28
    background::AUTO = 2e18
end
@xport @AutoParm struct CellDensity <: Density
    val::AUTO = 1e28
end
@xport @AutoParm struct ZeroField <: Field
    val::AUTO = 0.
end
@xport @AutoParm struct UniformField <: Field
    val::AUTO = 1.
end
@xport @AutoParm struct UniformTemp <: Temp
    val::AUTO = 0.
end

@xport @AutoParm struct UniformIonisation <: Ionisation end

@xport @AutoParm struct RatioIonisation <: Ionisation
    ratio::AUTO = (e) -> 0.5
end

@xport @AutoParm struct FunctionalIonisation <: Ionisation
    func::AUTO = (e) -> 0.5
end

@xport function FunctionalIonisationInit(;func=(e,E2) -> 1.,thresh=10.)
    cbins = [0.,(10. .^(-3:0.1:6.))...] #TODO: This is a dubplicate and depends on the initParams funciton
    cbins2 = [0.,(10. .^(-3:0.1:6.))...]
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
    return FunctionalIonisation(spl3)#, matrix,matrix2
end


@xport @AutoParm struct DiscreteIonisation <: Ionisation
    loss::AUTO = 5.
end

@xport @AutoParm struct DiscreteFunctionalIonisation <: Ionisation
    loss::AUTO = (e) -> 5.
end

@xport struct ZeroVelocity <: Velocity end

@xport @AutoParm struct InitVelocity <: Velocity
    ε::AUTO = 100.
    axis::AUTO = [0,0,1]
end
@xport @AutoParm struct GaussianEnergyVelocity <: Velocity
    ε::AUTO = 10.
    Δε_source::AUTO = 0.5
    axis::AUTO = [0,0,1]
end

@xport @AutoParm struct GaussianDirectionVelocity <: Velocity
    ε::AUTO = 10.
end

"""
MCMODEL

MCMODEL is the main constructor for MCSim.

```julia
    type = SWARM()                   # Used to set up specific models: LμJ(), SWARM() or CELL().
    particle_type = :electron        # Change particle type: :electron or :positron. This will modify ionisation behaviour.
    n = 100                          # Number of particles to simulate
    v_0 = ZeroVelocity()             # Initial velocity: ZeroVelocity()...
    r_0 = ZeroSpacial()              # Initial spacial: ZeroSpacial()...
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
    r_0::AUTO = ZeroSpacial()
    n0::AUTO = UniformDensity()
    En0::AUTO = UniformField()
    T::AUTO = UniformTemp()
    ion_type::AUTO = UniformIonisation()
    meas::AUTO = [:t]
    optimise_spatial::AUTO = true
end

@xport @AutoParm mutable struct SWARM
    t_f::AUTO = 2e-5t
end

"""
LμJ

This struct is used to initilize a micro-jet simulation. Used in `MCMODEL(type=LμJ())`.

Parameters:
```julia
    Δε_detector = 0.9       # Energy bin width for the detector
    Δχ_detector = 2         # Angle bin width for the detector in degrees
    χ = Inf                 # Only record measurements for a specific angle += resolution
    thresh = 1.             # Set the minimum energy to simulate electrons too
    jet_width = 15e-6       # Diameter for the liquid jet
    jet_center = [0.,0.,0.] # Set's an offset for the jet_center 
```

For more information, type ?MCMODEL
"""
@xport @AutoParm mutable struct LμJ
    Δε_detector::AUTO = 0.9
    Δχ_detector::AUTO = 2
    χ::AUTO = Inf
    thresh::AUTO = 1.
    jet_width::AUTO = 15e-6
    jet_center::AUTO = [0.,0.,0.]
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
    # Defined spacial parameters for the CELL
    Δε_detector::AUTO = 0.9
    thresh::AUTO = 1.
    z::AUTO = 15e-6
    z_n::AUTO = 100
    x::AUTO = 15e-6
    y::AUTO = 15e-6
    bin_no::AUTO = 100
    r::AUTO = 1e-6
end

@xport function LμJ(;jet_center=[0.,0.,0.],Δε_detector=0.9,Δχ_detector=2,χ=Inf,thresh=1.,jet_width=15e-6)
    LμJ(Δε_detector,Δχ_detector,χ,thresh,jet_width,jet_center)
end

@xport function CELL(;Δε_detector=0.9,thresh=1.,z=15e-6,z_n=100,x=15e-6,y=15e-6,bin_no=100,r=1e-6)
    CELL(Δε_detector,thresh,z,z_n,x,y,bin_no,r)
end

@AutoParm mutable struct MCPARAMS
    cbins::AUTO = [0.,(10. .^(-3:0.01:2))...]
    ebins::AUTO = nothing
    gas::AUTO = nothing
    splat_gas::AUTO = nothing
    total_freq::AUTO = nothing
    a::AUTO = nothing
    # n0::Union{Vector{Float64},Float64,Nothing} = nothing
    n0::AUTO = nothing
    rbins::Union{Vector{Float64},Nothing} = nothing
    χbins::Union{Vector{Float64},Nothing} = nothing
    tbins::Union{Vector{Float64},Nothing} = collect(LinRange(0,1e-5,10))
    zbins::Union{Vector{Float64},Nothing} = nothing
    xbins::Union{Vector{Float64},Nothing} = nothing
    ybins::Union{Vector{Float64},Nothing} = nothing
    length_scale::AUTO = nothing
    time_scale::AUTO = nothing
    energy_scale::AUTO = nothing

end

function MCPARAMS(;cbins=nothing,ebins=nothing,gas=nothing,splat_gas=nothing,total_freq=nothing,a=nothing,n0=nothing,rbins=nothing,χbins=nothing,tbins=nothing,zbins=nothing,xbins=nothing,ybins=nothing,length_scale=noting,time_scale=nothing,energy_scale=nothing)
    MCPARAMS(cbins,ebins,gas,splat_gas,total_freq,a,n0,rbins,χbins,tbins,zbins,xbins,ybins,length_scale,time_scale,energy_scale)
end


@xport @AutoParm mutable struct PART
    type::Symbol = :electron
    alive::Bool = true
    v::SArray{Tuple{3},Float64,1,3} = [0.,0.,0.]
    r::SArray{Tuple{3},Float64,1,3} = [0.,0.,0.]
    t::Float64 = 0.
    weight::Float64 = 1.
    tbin::Int64 = 1
    cbin::Int64 = 1
    ebin::Int64 = 1
    zbin::Int64 = 1
    ybin::Int64 = 1
    xbin::Int64 = 1
    rbin::Int64 = 1
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
    MEAS(MEASLμJ(),MEASCELL(),MEASG(),MEAST(params),MEASS(),MEASC(params))
end

function MEAS(type::LμJ,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(model,params),MEASCELL(),MEASG(),MEAST(),MEASS(),MEASC())
end

function MEAS(type::CELL,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(),MEASCELL(model,params),MEASG(model,params),MEAST(),MEASS(),MEASC())
end

@xport @AutoParm mutable struct COLLCOUNT
    elastic::AUTO = 0
    excitation::AUTO = 0
    ionisation::AUTO = 0
    attachment::AUTO = 0
end

@xport function sumColl(collcount::COLLCOUNT)
    return collcount.elastic+collcount.excitation+collcount.ionisation+collcount.attachment
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
    
    if :εχ ∈ model.meas
        εχ = zeros(length(params.ebins),length(params.χbins))
        εχ_reject = zeros(length(params.ebins),length(params.χbins))
    end

    return MEASLμJ(ε, χ, εχ,εχ_reject,zyt)

end

@xport @AutoParm mutable struct MEASCELL
    ε::AUTO = nothing
    εz::AUTO = nothing
    x::AUTO = nothing
    y::AUTO = nothing
    z::AUTO = nothing
    xy::AUTO = nothing
    xybins::AUTO = nothing
end

function MEASCELL(model::MCMODEL,params::MCPARAMS)
    ε, εz, x, y, z, xy, xybins = nothing, nothing, nothing, nothing, nothing, nothing, nothing
    if :ε ∈ model.meas
        ε = zeros(length(params.ebins))
    end
    if :εz ∈ model.meas
        εz = zeros(length(params.ebins),length(params.zbins))
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
    return MEASCELL(ε,εz,x,y,z,xy,xybins)
end

@AutoParm mutable struct MEASG
    coschi::AUTO = 0.
    coll_count::AUTO = COLLCOUNT()
    SP_ε::AUTO = nothing
    SP_r::AUTO = nothing
    SP_εs::AUTO = nothing
end

function MEASG(model::MCMODEL,params::MCPARAMS)
    coschi, coll_count, SP_ε, SP_r, SP_εs, SP_zs = 0., COLLCOUNT(), nothing, nothing, nothing, nothing
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

    MEAST(zeros(length(params.tbins)-1),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1),
          zeros(length(params.tbins)-1)
          )

end

@AutoParm mutable struct MEASS

end

@AutoParm mutable struct MEASC
    εᵢ::AUTO = nothing
    r::AUTO = nothing
    v::AUTO = nothing
    rv::AUTO = nothing
    weight::AUTO = nothing
end

function MEASC(params)
    return MEASC(zeros(length(params.tbins)-1),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1)
          )
end


@AutoParm @xport mutable struct SWARMMEAS
    ε::Float64 = 0.
    εs::Vector{Float64} = [0.]
    bW::Float64 = 0.
    bWs::Vector{Float64} = [0.]
    fW::Float64 = 0.
    fWs::Vector{Float64} = [0.]
    bDt::Float64 = 0.
    bDts::Vector{Float64} = [0.]
    fDt::Float64 = 0.
    fDts::Vector{Float64} = [0.]
    bDl::Float64 = 0.
    bDls::Vector{Float64} = [0.]
    fDl::Float64 = 0.
    fDls::Vector{Float64} = [0.]
end

@fastmath function findMinInd(x::Array{Float64,1},y::Float64)::Int64
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
