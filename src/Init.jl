module Init

using ..Structs
using PhysicalConstants.CODATA2018, Unitful, LinearAlgebra, StaticArrays, AndExport, Gases, LXCatParser, AutoParameters, Interpolations, QuadGK

export initParams, initMeas, PART, getEnergy

const amu = ustrip(float(AtomicMassConstant))
const q = ustrip(float(ElementaryCharge))
const me = ustrip(float(ElectronMass))
const k = ustrip(float(BoltzmannConstant))
const ħ = ustrip(float(ReducedPlanckConstant))
const Td = 1e-21

@AutoParm struct MCGAS
    cs::AUTO
    dcs::AUTO
    type::AUTO
    thresh::AUTO
    m::AUTO 
    w::AUTO
end


function initParams(input_gas=MaxwellModel(),model::MCMODEL=MCMODEL())

    cbins = [0.,(10. .^(-3:0.1:2.))...]
    tbins, ebins, rbins, χbins, zbins, ybins, xbins = generateBins(model)
    n0=buildDensity(model.n0,tbins)
    a=buildField(model.En0)

    length_scale = model.n0.val*1e-20
    time_scale = model.n0.val*1e-20*sqrt(2*q/me)
    energy_scale = (length_scale/time_scale)^2

    splat_gas=splatCS(input_gas,energy_scale,cbins)
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

@xport function unscale!(type::SWARM,params,meas::MEAS)
    ls = 1/params.length_scale
    ts = 1/params.time_scale
    es = 1/params.energy_scale
    meas.meast.εᵢ .*= es
    meas.meast.v² .*= es
    meas.meast.r² .*= ls^2
    meas.meast.r .*= ls
    meas.meast.v .*= ls/ts
    meas.meast.rv .*= ls^2/ts
end

@xport function unscale!(type::LμJ,params,meas::MEAS)
    ls = 1/params.length_scale
    ts = 1/params.time_scale
    es = 1/params.energy_scale
end

function initMeas(model,params::MCPARAMS = MCPARAMS())
    MEAS(model.type,model,params)
end

#### Frequency ####

function buildFreq(model,gas,cbins,energy_scale)

    total_freq = zeros(length(cbins))
    en_spacing = cbins[end]-cbins[end-1]
    for i in 1:length(total_freq)-1
        total_freq[i] = maximum([sum([gas.cs[i](en)*sqrt(2*en*q/me)*gas.w[i] for i in eachindex(gas.cs)]) for en in LinRange(cbins[i]+1e-6,cbins[i+1],1000)])*1e-20
        # total_freq[i] = maximum([sum([cs(en) for (cs,type,thresh) in gas]) for en in LinRange(cbins[i]+0.000001,cbins[i+1],1000)])*1e-20
        #NOTE: This is no longer a frequency
    end
    total_freq[end] = maximum([sum([gas.cs[i](en)*sqrt(2*en*q/me)*gas.w[i] for i in eachindex(gas.cs)]) for en in LinRange(cbins[end],cbins[end]+en_spacing,1000)])*1e-20
    return total_freq
end

@xport function splatCS(gas,energy_scale,cbins)
    gases = [gas.elastic...,gas.excitation...,gas.ionisation...,gas.attachment...]
    gases_dcss = [gas.elastic_DCS...,gas.excitation_DCS...,gas.ionisation_DCS...,gas.attachment_DCS...] 

    en_spacing = cbins[end]-cbins[end-1]
    ens = LinRange(1e-6,cbins[end]+en_spacing+1.,100)
    ens2 = 1e-6:1e-4:(cbins[end]+en_spacing+1.)
    θs = LinRange(1e-6,π,100)

    gases_itp::Vector{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Linear{Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}}, BSpline{Linear{Throw{OnGrid}}}, Flat{Nothing}}} = []
    for gas in gases
        push!(gases_itp, LinearInterpolation(ens2,gas.(ens2),extrapolation_bc=Flat()))
    end

    dcss::Vector{Vector{Interpolations.Extrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Float64, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Tuple{Vector{Float64}}}, Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}}, Interpolations.Flat{Nothing}}}} = []
    # replace_nan(v) = map(x -> isnan(x) ? zero(x)+1e-10 : x, v)
    # replace_0(v) = map(x -> x == 0 ? zero(x)+1e-10 : x, v)
    for σ in gases_dcss
        DCS_matrix::Matrix{Float64} = zeros(length(ens),length(θs))
        for i in eachindex(ens)
            func(chi) = σ(ens[i],chi)*sin(chi) < 1e-10 ? 1e-10 : σ(ens[i],chi)*sin(chi)
            norm = quadgk(func,0,pi)[1]
            for j in eachindex(θs)
                if j == 1
                    DCS_matrix[i,j] = quadgk(func,0,θs[j])[1]/norm
                elseif j == length(θs)
                    DCS_matrix[i,j] = 1. - sum(DCS_matrix[i,:])
                else
                    DCS_matrix[i,j] = quadgk(func,θs[j-1],θs[j])[1]/norm
                end
            end
        end
        # DCS_matrix = replace_0(replace_nan(DCS_matrix))
        DCS_matrix = cumsum(DCS_matrix,dims=2)
        spl2 = LinearInterpolation((ens,θs),DCS_matrix,extrapolation_bc=Flat())	#2D spline of the matrix
        # func3(e,R) = LinearInterpolation(spl2(e,θs),cos.(θs))(R)::Float64	#Interpolates the above vector and evaluates a value of chi for a given R
        # @show typeof(LinearInterpolation(spl2(1.,θs),cos.(θs)))	#Interpolates the above vector and evaluates a value of chi for a given R
        #TODO: Try this with energy dependency aswell
        push!(dcss,[LinearInterpolation(Interpolations.deduplicate_knots!(spl2(c,θs)),cos.(θs),extrapolation_bc=Flat()) for c in cbins])
    end
    
    
    types = [[:ELASTIC for i in 1:length(gas.elastic)]...,[:EXCITATION for i in 1:length(gas.excitation)]...,[:IONIZATION for i in 1:length(gas.ionisation)]...,[:ATTACHMENT for i in 1:length(gas.attachment)]...]

    threshes = [[0. for thresh in gas.elastic]...,[thresh for thresh in gas.exlosses]...,[thresh for thresh in gas.iolosses]...,[0. for thresh in gas.attachment]...].*energy_scale .*q

    masses = [[m for m in gas.Mₑₗ]...,[m for m in gas.Mₑₓ]...,[m for m in gas.Mₐₜ]...,[m for m in gas.Mᵢₒ]...] .* amu

    ws = [[w for w in gas.wₑₗ]...,[w for w in gas.wₑₓ]...,[w for w in gas.wₐₜ]...,[w for w in gas.wᵢₒ]...]


    MCGAS(gases_itp,dcss,types,threshes,masses,ws)
    # return gases_itp,dcss,types,threshes,masses,ws
end

#### * Bins
# This tries a best guess at what parameters will be needed but is in no way good. Once params are initialized you will need to customize the ones you want as you need

function generateBins(model)

    m = string.(model.meas)

    tbins = 1 ∈ occursin.("t",m) ? buildTBins(model.type,model) : nothing
    ebins = 1 ∈ occursin.("e",m) ? buildEnBins(model.type,model) : nothing
    zbins = 1 ∈ occursin.("z",m) ? buildZBins(model.type,model) : nothing
    ybins = 1 ∈ occursin.("y",m) ? buildYBins(model.type,model) : nothing
    xbins = 1 ∈ occursin.("x",m) ? buildXBins(model.type,model) : nothing
    rbins = 1 ∈ occursin.("r",m) ? buildRBins(model.type,model,model.n0) : nothing
    χbins = 1 ∈ occursin.("χ",m) ? buildχBins(model.type,model) : nothing

    return tbins, ebins, rbins, χbins, zbins, ybins, xbins

end

function buildRBins(type::SWARM,model,n0::UniformDensity)
    return LinRange(0,1e-6,10)
end

function buildRBins(type::LμJ,model,n0::LμJDensity)
    return [0.,(type.jet_width/2)^2,1e20]
end

function buildEnBins(type::SWARM,model)
    return 10 .^LinRange(-3,2,1000)
end

function buildEnBins(type::LμJ,model)
    return collect(model.type.thresh:model.type.Δε_detector:(model.v_0.ε+2*model.v_0.Δε_source))
end

function buildTBins(type::LμJ,model)
    return LinRange(3e-13,320e-13,10)
end

function buildTBins(type::SWARM,model)
    return LinRange(0,type.t_f,1000)
end

function buildχBins(type::LμJ,model)
    return LinRange(0,π,100)
end

function buildχBins(type::SWARM,model)
    return LinRange(0,π,100)
end

function buildZBins(type::LμJ,model)
    return LinRange(-7.5e-6,-6.5e-6,200)
end

function buildZBins(type::SWARM,model)
    return LinRange(-7.5e-6,-6.5e-6,200)
end

function buildYBins(type::LμJ,model)
    return LinRange(-7e-7,7e-7,200)
end

function buildYBins(type::SWARM,model)
    return LinRange(-7e-7,7e-7,200)
end

function buildXBins(type::LμJ,model)
    return LinRange(-7e-7,7e-7,200)
end

function buildXBins(type::SWARM,model)
    return LinRange(-7e-7,7e-7,200)
end

#### Density ####

function buildDensity(n0::UniformDensity,rbins)
    n0.val
end

function buildDensity(n0::LμJDensity,rbins)
    [n0.val,0.,0.]
end
#### Electric Field ####

function buildField(En0::UniformField)
    SVector(0.,0.,q*Td*En0.val/me)
end



#### * PARTICLE

function Structs.PART(type::SWARM,model,params)

    particle = PART()
    particle.r = initPartr(model.r_0,params).*params.length_scale
    particle.v = initPartv(model.v_0,params).*params.length_scale ./params.time_scale
    particle.cbin = findMinInd(params.cbins,getEnergy(particle.v))
    if params.zbins != nothing
        particle.zbin = findMinInd(params.zbins,particle.r[3])
    end
    
    if params.ybins != nothing
        particle.ybin = findMinInd(params.ybins,particle.r[2])
    end
        

    return particle

end
function Structs.PART(type::LμJ,model,params)

    particle = PART()
    particle.r = initPartr(model.r_0,params).*params.length_scale
    particle.v = initPartv(model.v_0,params).*params.length_scale ./params.time_scale
    particle.cbin = findMinInd(params.cbins,getEnergy(particle.v))
    particle.rbin = findMinInd(params.rbins,dot(particle.r,particle.r))

    if params.zbins != nothing
        particle.zbin = findMinInd(params.zbins,particle.r[3])
    end
    
    if params.ybins != nothing
        particle.ybin = findMinInd(params.ybins,particle.r[2])
    end
    return particle

end

function getEnergy(v)
    0.5*me*dot(v,v)/q
end

function getVelocity(ε)
    sqrt(ε*2*q/me)
end


function initPartr(r_0::ZeroSpacial,params)
    [0.,0.,0.]
end
function initPartr(r_0::InitSpacial,params)
    r_0.axis .* r_0.val
end
function initPartr(r_0::GaussianLμJSpacial,params)
    SD = r_0.Δy_source/(2*sqrt(2*log(2)))

    val = 0.

    R=randn()
    # while val < -SD*R
        # R=randn()
    # end
    val .+ SD*R .* r_0.gaussian_axis .+ r_0.val .* r_0.axis
end
function initPartv(v_0::ZeroVelocity,params)
    [0.,0.,0.]
end

function initPartv(v_0::InitVelocity,params)
    v_0.axis .* getVelocity(v_0.ε)
end

function initPartv(v_0::GaussianVelocity,params)
    SD = v_0.Δε_source/(2*sqrt(2*log(2)))
    electron_v2eV = 0.5 * me / q

    ea = v_0.ε

    R=randn()
    while ea < -SD*R
        R=randn()
    end
    new_ea=ea+SD*R

    sqrt.(new_ea*v_0.axis/electron_v2eV)
end


end
