module Init

using ..Utils
using ..Structs
using LinearAlgebra, StaticArrays, AndExport, Gases, LXCatParser, AutoParameters, Interpolations, QuadGK, MCSimUnits

export initParams, initMeas, PART, getEnergy, test

@AutoParm struct MCGAS
    cs::AUTO
    dcs::AUTO
    type::AUTO
    thresh::AUTO
    m::AUTO
    w::AUTO
end


function initParams(input_gas=MaxwellModel(),model::MCMODEL=MCMODEL())

    cbins = [0.,(10. .^(-6:0.1:6.))...] .*eV
    # cbins = collect(LinRange(0,101,100))
    tbins, ebins, rbins, χbins, zbins, ybins, xbins = generateBins(model)
    n0=buildDensity(model.n0,rbins,zbins)
    a=buildField(model.En0)

    length_scale = model.n0.val*1e-20
    time_scale = model.n0.val*1e-20*sqrt(2*q/me)
    energy_scale = (length_scale/time_scale)^2

    splat_gas=splatCS(input_gas,model.T.val,energy_scale,ustrip.(cbins))
    total_freq=buildFreq(model,splat_gas,cbins,energy_scale)

    params = MCPARAMS(
        ebins=ebins,
        tbins=tbins,
        zbins=zbins,
        ybins=ybins,
        xbins=xbins,
        cbins=cbins,
        rbins=rbins,
        χbins=χbins,
        gas=input_gas,
        splat_gas=splat_gas,
        n0=n0,
        a=a,
        total_freq=total_freq,
        length_scale = length_scale,
        time_scale = time_scale,
        energy_scale = energy_scale
        )

    return params
end

# * Scale functions
@xport function scale!(model,type::LμJ,params)
    ls = params.length_scale
    ts = params.time_scale
    es = params.energy_scale

    model.type.jet_center .*= ls
end

@xport function scale!(model,type::CELL,params)
    ls = params.length_scale
    ts = params.time_scale
    es = params.energy_scale
end

@xport function scale!(model,type::SWARM,params)
    ls = params.length_scale
    ts = params.time_scale
    es = params.energy_scale
end
@xport function scale!(params)
    ls = params.length_scale
    ts = params.time_scale
    es = params.energy_scale

    params.cbins *= es
    params.ebins *= es
    params.total_freq /= ts#/ls^3
    params.a *= ls^4/ts^2
    params.n0 /= ls^3
    params.rbins *= ls^2
    params.xbins *= ls
    params.ybins *= ls
    params.zbins *= ls
    params.tbins *= ts

    return params
end

@xport function unscale!(params)
    ls = 1/params.length_scale
    ts = 1/params.time_scale
    es = 1/params.energy_scale

    params.cbins *= es
    params.ebins *= es
    params.total_freq /= ts#/ls^3
    params.a *= ls^4/ts^2
    params.n0 /= ls^3
    params.rbins *= ls^2
    params.xbins *= ls
    params.ybins *= ls
    params.zbins *= ls
    params.tbins *= ts

    return params
end

@xport function unscale!(type::SWARM,model,params,meas::MEAS)
    ls = 1/params.length_scale
    ts = 1/params.time_scale
    es = 1/params.energy_scale
    meas.meast.εᵢ .*= es
    meas.meast.v² .*= es
    meas.meast.r² .*= ls^2
    meas.meast.r .*= ls
    meas.meast.v .*= ls/ts
    meas.meast.rv .*= ls^2/ts
    meas.measc.r .*= ls
    meas.measc.v .*= ls/ts
    meas.measc.rv .*= ls^2/ts
end

@xport function unscale!(type::LμJ,model,params,meas::MEAS)
    ls = 1/params.length_scale
    ts = 1/params.time_scale
    es = 1/params.energy_scale
    model.type.jet_center .*= ls
end

@xport function unscale!(type::CELL,model,params,meas::MEAS)
    ls = 1/params.length_scale
    ts = 1/params.time_scale
    es = 1/params.energy_scale
    meas.measg.SP_ε *= es
    :SP_z in model.meas ? meas.measg.SP_εs .*= es : nothing
    meas.measg.SP_r *= ls
end

function initMeas(model,params::MCPARAMS = MCPARAMS())
    MEAS(model.type,model,params)
end

#### Frequency ####

function buildFreq(model,gas,cbins,energy_scale)
    total_freq = zeros(length(cbins)) .*m^3 ./s
    en_spacing = cbins[end]-cbins[end-1]
    for i in 1:length(total_freq)-1
        total_freq[i] = maximum([sum([gas.cs[i](en)*sqrt(2*en/me)*gas.w[i] for i in eachindex(gas.cs)]) for en in LinRange(cbins[i]+1e-6eV,cbins[i+1],1000)])
        # total_freq[i] = maximum([sum([cs(en) for (cs,type,thresh) in gas]) for en in LinRange(cbins[i]+0.000001,cbins[i+1],1000)])*1e-20
        #NOTE: This is no longer a frequency
    end
    total_freq[end] = maximum([sum([gas.cs[i](en)*sqrt(2*en/me)*gas.w[i] for i in eachindex(gas.cs)]) for en in LinRange(cbins[end],cbins[end]+en_spacing,1000)])
    return total_freq
end

@xport function addSuperElastics(gas,T)
    excitations_super::Vector{Function} = []
    excitations::Vector{Function} = []

    # Z_total = Z_T(gas.exlosses,ustrip(T))

    no_gases = maximum(gas.ex_nos)
    for j in 1:no_gases
        Z_total = 1.
        for (i,cs) in enumerate(gas.excitation) if gas.excitation[i].type === :ROTATION && gas.ex_nos[i] == j
            Z_total += Z(gas.exlosses[i],ustrip(T))
        end end
        for (i,cs) in enumerate(gas.excitation) if gas.excitation[i].type === :ROTATION && gas.ex_nos[i] == j
            push!(excitations_super, (ε)->(ε + gas.exlosses[i])*cs(ε+gas.exlosses[i])*Z(gas.exlosses[i],ustrip(T))./(Z_total * ε))
            push!(excitations, (ε)->cs(ε)./(Z_total))
        end end

        for t in (:VIBRATION, :EXCITATION, :DISSOCIATION)
            for (i,cs) in enumerate(gas.excitation) if gas.excitation[i].type === t && gas.ex_nos[i] == j
                push!(excitations_super, (ε)->(ε + gas.exlosses[i])*cs(ε+gas.exlosses[i])*Z(gas.exlosses[i],ustrip(T))./((1+Z(gas.exlosses[i],ustrip(T))) * ε))
                push!(excitations, (ε)->cs(ε)./(1+Z(gas.exlosses[i],ustrip(T))))
            end end
        end
    end

    gases = [gas.elastic...,excitations...,excitations_super...,gas.ionisation...,gas.attachment...]

    gases_dcss = [gas.elastic_DCS...,gas.excitation_DCS...,gas.excitation_DCS...,gas.ionisation_DCS...,gas.attachment_DCS...]
    return gases, gases_dcss
end


@xport function splatCS(gas,T,energy_scale,cbins)
    gases, gases_dcss = addSuperElastics(gas,T)
    # gases = [gas.elastic...,gas.excitation...,gas.ionisation...,gas.attachment...]
    # gases_dcss = [gas.elastic_DCS...,gas.excitation_DCS...,gas.ionisation_DCS...,gas.attachment_DCS...]
    en_spacing = cbins[end]-cbins[end-1]
    # ens = LinRange(1e-6,cbins[end]+en_spacing+1.,100)
    ens = cbins[2:end]
    # ens2 = 1e-6:1e-4:(cbins[end]+en_spacing+1.)
    ens2 = cbins[2:end]
    θs = LinRange(1e-6,π,1000)

    gases_itp = [LinearInterpolation(ens2.*eV,gas.(ens2).*u"Å^2",extrapolation_bc=Flat()) for gas in gases]

    dcss::Vector{Vector{Interpolations.Extrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Vector{Float64}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Interpolations.Flat{Nothing}}}} = []
    # replace_nan(v) = map(x -> isnan(x) ? zero(x)+1e-10 : x, v)
    # replace_0(v) = map(x -> x == 0 ? zero(x)+1e-10 : x, v)
    for σ in gases_dcss
        DCS_matrix::Matrix{Float64} = zeros(length(ens),length(θs))
        for i in eachindex(ens)
            func(chi)::Float64 = σ(ens[i],chi)*sin(chi) < 1e-10 ? 1e-10 : σ(ens[i],chi)*sin(chi)
            norm = quadgk(func,0,pi)[1]::Float64
            for j in eachindex(θs)
                if j == 1
                    DCS_matrix[i,j] = quadgk(func,0,θs[j])[1]::Float64/norm
                elseif j == length(θs)
                    DCS_matrix[i,j] = 1. - sum(DCS_matrix[i,:])
                else
                    DCS_matrix[i,j] = quadgk(func,θs[j-1],θs[j])[1]::Float64/norm
                end
            end
        end
        # DCS_matrix = replace_0(replace_nan(DCS_matrix))
        DCS_matrix = cumsum(DCS_matrix,dims=2)
        spl2 = LinearInterpolation((ens,θs),DCS_matrix,extrapolation_bc=Flat())	#2D spline of the matrix
        # func3(e,R) = LinearInterpolation(spl2(e,θs),cos.(θs))(R)::Float64	#Interpolates the above vector and evaluates a value of chi for a given R
        # @show typeof(LinearInterpolation(spl2(1.,θs),cos.(θs)))	#Interpolates the above vector and evaluates a value of chi for a given R
        #TODO: Try this with energy dependency aswell
        push!(dcss,[LinearInterpolation(Interpolations.deduplicate_knots!(spl2(c,θs),move_knots=true),cos.(θs),extrapolation_bc=Flat()) for c in cbins])
    end


    types = [[:ELASTIC for i in 1:length(gas.elastic)]...,[:EXCITATION for i in 1:length(gas.excitation)]...,[:EXCITATION for i in 1:length(gas.excitation)]...,[:IONIZATION for i in 1:length(gas.ionisation)]...,[:ATTACHMENT for i in 1:length(gas.attachment)]...]

    threshes = [[0. for thresh in gas.elastic]...,[thresh for thresh in gas.exlosses]...,[-thresh for thresh in gas.exlosses]...,[thresh for thresh in gas.iolosses]...,[0. for thresh in gas.attachment]...].*eV

    masses = [[m for m in gas.Mₑₗ]...,[m for m in gas.Mₑₓ]...,[m for m in gas.Mₑₓ]...,[m for m in gas.Mₐₜ]...,[m for m in gas.Mᵢₒ]...] .* amu

    ws = [[w for w in gas.wₑₗ]...,[w for w in gas.wₑₓ]...,[w for w in gas.wₑₓ]...,[w for w in gas.wₐₜ]...,[w for w in gas.wᵢₒ]...]


    MCGAS(gases_itp,dcss,types,threshes,masses,ws)
    # return gases_itp,dcss,types,threshes,masses,ws
end

#### * Bins
# This tries a best guess at what parameters will be needed but is in no way good. Once params are initialized you will need to customize the ones you want as you need

function generateBins(model)

    m = string.(model.meas)

    tbins = 1 ∈ occursin.("t",m) ? buildTBins(model.type,model) : nothing
    ebins = 1 ∈ occursin.("ε",m) ? buildEnBins(model.type,model) : nothing
    zbins = 1 ∈ occursin.("z",m) ? buildZBins(model.type,model) : nothing
    ybins = 1 ∈ occursin.("y",m) ? buildYBins(model.type,model) : nothing
    xbins = 1 ∈ occursin.("x",m) ? buildXBins(model.type,model) : nothing
    rbins = 1 ∈ occursin.("r",m) ? buildRBins(model.type,model,model.n0) : nothing
    χbins = 1 ∈ occursin.("χ",m) ? buildχBins(model.type,model) : nothing

    return tbins, ebins, rbins, χbins, zbins, ybins, xbins

end

# TODO: Combine these into one function per model

function buildRBins(type::SWARM,model,n0::UniformDensity)
    return LinRange(0m,1e-6m,10)
end

function buildEnBins(type::SWARM,model)
    return 10 .^LinRange(-3,2,1000) .*eV
end

function buildTBins(type::SWARM,model)
    return LinRange(type.t_f/2,type.t_f,100)
end

function buildχBins(type::SWARM,model)
    return LinRange(0,π,90) .* u"rad"
end

function buildZBins(type::SWARM,model)
    return LinRange(-7.5e-6,-6.5e-6,200) .*m
end

function buildYBins(type::SWARM,model)
    return LinRange(-7e-7,7e-7,200) .*m
end

function buildXBins(type::SWARM,model)
    return LinRange(-7e-7,7e-7,200) .*m
end

function buildRBins(type::CELL,model,n0)
    return LinRange(0,type.r,type.bin_no) .*m^2
end

function buildEnBins(type::CELL,model)
    #TODO: Make specific to each velocity type
    #return collect(model.type.thresh:model.type.Δε_detector:(model.v_0.ε+2*model.v_0.Δε_source))
    return collect(model.type.thresh:10:model.v_0.ε)
end

function buildTBins(type::CELL,model)
    return LinRange(0s,type.t_f,1000)
end

function buildχBins(type::CELL,model)
    return LinRange(0,π,100) .*u"rad"
end

function buildZBins(type::CELL,model)
    return LinRange(0. *m,type.z,type.z_n)
end

function buildYBins(type::CELL,model)
    return LinRange(-type.y,type.y,200)
end

function buildXBins(type::CELL,model)
    return LinRange(-type.x,type.x,200)
end

function buildRBins(type::LμJ,model,n0::LμJDensity)
    return [0. *m^2,(type.jet_width/2)^2,0.02^2 *m^2]
end

function buildRBins(type::LμJ,model,n0::LμJDensityGradient)
    return adaptiveMesh(n0.func,0. *m,0.02 *m,n0.δtol).^2
end

function buildEnBins(type::LμJ,model)
    return collect(model.type.thresh:model.type.Δε_detector:(model.v_0.ε+2*model.v_0.Δε_source))
end

function buildTBins(type::LμJ,model)
    return LinRange(3e-13,320e-13,10) .*s
end

function buildχBins(type::LμJ,model)
    return LinRange(0,π,100) .*u"rad"
end

function buildZBins(type::LμJ,model)
    return LinRange(-7.5e-6,-6.5e-6,200) .*m
end

function buildYBins(type::LμJ,model)
    return LinRange(-7e-7,7e-7,200) .*m
end

function buildXBins(type::LμJ,model)
    return LinRange(-7e-7,7e-7,200) .*m
end

#### Density ####

function buildDensity(n0::UniformDensity,rbins,zbins)
    n0.val
end

function buildDensity(n0::LμJDensityGradient,rbins,zbins)
    n0.func.(sqrt.(rbins))
end

function buildDensity(n0::LμJDensity,rbins,zbins)
    [n0.val,n0.background,0. *m^-3]
end

function buildDensity(n0::CellDensity,rbins,zbins)
    n0.val .* ones(length(zbins))
end

#### Electric Field ####

function buildField(En0::UniformField)
    SVector(0. *m^4*s^-2,0. *m^4*s^-2,upreferred((q*En0.val/me)))
end

function buildField(En0::ZeroField)
    SVector(0. *m^4/s^-2,0. *m^4/s^-2,0. *m^4/s^-2)
end


#### * PARTICLE

function Structs.PART(model,params)

    particle = PART()
    particle.type = model.particle_type
    particle.r = initPartr(model.r_0,params)#.*params.length_scale
    particle.v = initPartv(model.v_0,params)#.*params.length_scale ./params.time_scale
    particle.cbin = findMinInd(params.cbins,getEnergy(particle.v))
    if params.ebins != nothing
        particle.ebin = findMinInd(params.ebins,getEnergy(particle.v))
    end

    if params.rbins != nothing
        particle.rbin = findMinInd(params.rbins,dot(particle.r,particle.r))
    end

    if params.zbins != nothing
        particle.zbin = findMinInd(params.zbins,particle.r[3])
    end

    if params.xbins != nothing
        particle.xbin = findMinInd(params.xbins,particle.r[1])
    end

    if params.ybins != nothing
        particle.ybin = findMinInd(params.ybins,particle.r[2])
    end


    return particle

end



function initPartr(r_0::ZeroSpatial,params)
    r = SArray{Tuple{3},Float64,1,3}(0.,0.,0.)
    return r*m
end
function initPartr(r_0::InitSpatial,params)
    r = r_0.axis .* r_0.val
end
function initPartr(r_0::GaussianLμJSpatial,params)
    SD = r_0.Δy_source/(2*sqrt(2*log(2)))

    val = 0.0m

    R=randn()
    # while val < -SD*R
        # R=randn()
    # end
    val .+ SD*R .* r_0.gaussian_axis .+ r_0.val .* r_0.axis
end
function initPartv(v_0::ZeroVelocity,params)
    r = SArray{Tuple{3},Float64,1,3}(0.,0.,0.)
    return r*m/s
end

function initPartv(v_0::InitVelocity,params)
    #TODO: Make all these static arrays
    v= v_0.axis .* getVelocity(v_0.ε)
    return v
end

function initPartv(v_0::GaussianDirectionVelocity,params)
    electron_v2eV = 0.5 * me
    coschi = 1-2*rand()
    sinchi = sqrt(1-coschi^2)
    phi = rand()*2*pi

    sqrt(v_0.ε/electron_v2eV) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
end


function initPartv(v_0::GaussianEnergyVelocity,params)
    SD = v_0.Δε_source/(2*sqrt(2*log(2)))
    electron_v2eV = 0.5 * me

    ea = v_0.ε

    R=randn()
    while ea < -SD*R
        R=randn()
    end
    new_ea=ea+SD*R

    sqrt.(new_ea/electron_v2eV)*v_0.axis
end


end
