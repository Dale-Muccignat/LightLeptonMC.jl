module Init

using ..Structs
using PhysicalConstants.CODATA2018, Unitful, LinearAlgebra, StaticArrays, AndExport, Gases, LXCatParser, AutoParameters

export initParams, initMeas, PART, getEnergy

const amu = ustrip(float(AtomicMassConstant))
const q = ustrip(float(ElementaryCharge))
const me = ustrip(float(ElectronMass))
const k = ustrip(float(BoltzmannConstant))
const ħ = ustrip(float(ReducedPlanckConstant))
const Td = 1e-21

@AutoParm struct MCGAS
    cs::AUTO
    type::AUTO
    thresh::AUTO
    m::AUTO 
    w::AUTO
end


function initParams(type::SWARM,input_gas=MaxwellModel,model::MCMODEL=MCMODEL())

    ebins = [0.,(10. .^(-3:0.1:2.))...]
    tbins = collect(LinRange(0,type.t_f,1000))
    n0=buildDensity(model.n0,tbins)
    a=buildField(model.En0)

    length_scale = model.n0.val*1e-20
    time_scale = model.n0.val*1e-20*sqrt(2*q/me)
    energy_scale = (length_scale/time_scale)^2
    splat_gas=splatCS(input_gas,energy_scale)
    total_freq=buildFreq(model,splat_gas,ebins,energy_scale)

    params = MCPARAMS(
        ebins=ebins,
        gas=input_gas,
        splat_gas=splat_gas,
        tbins=tbins,
        n0=n0,
        a=a,
        total_freq=total_freq,
        length_scale = length_scale,
        time_scale = time_scale,
        energy_scale = energy_scale
        )

    return params
end

@xport function scale!(params)
    ls = params.length_scale
    ts = params.time_scale
    es = params.energy_scale
    
    params.ebins *= es
    params.ebinsLμJ *= es
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
    
    params.ebins *= es
    params.ebinsLμJ *= es
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

function initParams(type::LμJ,input_gas,model::MCMODEL=MCMODEL())

    ebins = [0.,(10. .^(-3:0.1:log10(model.v_0.ε)+0.1))...]
    rbins = buildRBins(type,model.n0)
    n0=buildDensity(model.n0,rbins)
    a=buildField(model.En0)

    length_scale = model.n0.val*1e-20
    time_scale = model.n0.val*1e-20*sqrt(2*q/me)
    energy_scale = (length_scale/time_scale)^2

    splat_gas=splatCS(input_gas,energy_scale)
    total_freq=buildFreq(model,splat_gas,ebins,energy_scale)

    params = MCPARAMS(
        ebins=ebins,
        ebinsLμJ=collect(model.type.thresh:model.type.Δε_detector:(model.v_0.ε+5.)),
        gas=input_gas,
        splat_gas=splat_gas,
        rbins=rbins,
        n0=n0,
        a=a,
        total_freq=total_freq,
        length_scale = length_scale,
        time_scale = time_scale,
        energy_scale = energy_scale
        )

    return params
end
function initMeas(model,params::MCPARAMS = MCPARAMS())

    MEAS(model.type,model,params)

end

#### Frequency ####

function buildFreq(model,gas,ebins,energy_scale)

    total_freq = zeros(length(ebins))
    en_spacing = ebins[end]-ebins[end-1]
    for i in 1:length(total_freq)-1
        total_freq[i] = maximum([sum([gas.cs[i](en)*sqrt(2*en*q/me)*gas.w[i] for i in eachindex(gas.cs)]) for en in LinRange(ebins[i]+0.000001,ebins[i+1],1000)])*1e-20
        # total_freq[i] = maximum([sum([cs(en) for (cs,type,thresh) in gas]) for en in LinRange(ebins[i]+0.000001,ebins[i+1],1000)])*1e-20
        #NOTE: This is no longer a frequency
    end
    total_freq[end] = maximum([sum([gas.cs[i](en)*sqrt(2*en*q/me)*gas.w[i] for i in eachindex(gas.cs)]) for en in LinRange(ebins[end],ebins[end]+en_spacing,1000)])*1e-20
    return total_freq
end

@xport function splatCS(gas,energy_scale)
    gases = (gas.elastic...,gas.excitation...,gas.ionisation...,gas.attachment...)
    types = ([:ELASTIC for i in 1:length(gas.elastic)]...,[:EXCITATION for i in 1:length(gas.excitation)]...,[:IONIZATION for i in 1:length(gas.ionisation)]...,[:ATTACHMENT for i in 1:length(gas.attachment)]...)
    threshes = ([0. for thresh in gas.elastic]...,[thresh for thresh in gas.exlosses]...,[thresh for thresh in gas.iolosses]...,[0. for thresh in gas.attachment]...).*energy_scale .*q
    masses = ([m for m in gas.Mₑₗ]...,[m for m in gas.Mₑₓ]...,[m for m in gas.Mₐₜ]...,[m for m in gas.Mᵢₒ]...) .* amu
    ws = ([w for w in gas.wₑₗ]...,[w for w in gas.wₑₓ]...,[w for w in gas.wₐₜ]...,[w for w in gas.wᵢₒ]...)


    MCGAS(gases,types,threshes,masses,ws)
end

#### * RBins

function buildRBins(type::SWARM,n0::UniformDensity)
    return nothing
end

function buildRBins(type::LμJ,n0::LμJDensity)
    return [0.,(type.jet_width/2)^2,1e20]
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
    particle.ebin = findMinInd(params.ebins,getEnergy(particle.v))

    return particle

end
function Structs.PART(type::LμJ,model,params)

    particle = PART()
    particle.r = initPartr(model.r_0,params).*params.length_scale
    particle.v = initPartv(model.v_0,params).*params.length_scale ./params.time_scale
    particle.ebin = findMinInd(params.ebins,getEnergy(particle.v))
    particle.rbin = findMinInd(params.rbins,dot(particle.r,particle.r))

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

# @xport function getGas(input_gas,ebins)
#     CSs::Vector{CrossSection} = []
#     for gas in input_gas.elastic
#         push!(CSs,CrossSection(:ELASTIC,"",me/input_gas.Mₑₗ[1],0.,1.,ebins,gas.(ebins)))
#     end
#     for gas in input_gas.excitation
#         push!(CSs,CrossSection(:EXCITATION,"",me/input_gas.Mₑₗ[1],find_loss(gas),1.,ebins,gas.(ebins)))
#     end
#     for gas in input_gas.ionisation
#         push!(CSs,CrossSection(:IONIZATION,"",me/input_gas.Mₑₗ[1],find_loss(gas),1.,ebins,gas.(ebins)))
#     end
    
#     return CSs

# end

end
