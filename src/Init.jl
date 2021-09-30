module Init

using MCModels
using PhysicalConstants.CODATA2018, Unitful, LinearAlgebra, StaticArrays, AndExport

export initParams, initMeas, PART, getEnergy

const amu = ustrip(float(AtomicMassConstant))
const q = ustrip(float(ElementaryCharge))
const me = ustrip(float(ElectronMass))
const k = ustrip(float(BoltzmannConstant))
const ħ = ustrip(float(ReducedPlanckConstant))
const Td = 1e-21

function initParams(model::MCMODEL=MCMODEL())

    params = MCPARAMS()

    buildFreq!(model,params) #TODO
    buildDensity!(model.n0,params)
    buildField!(model.n0,model.En0,params)
    return params
end

function initMeas(params::MCPARAMS = MCPARAMS())

    MEAS(params)


end

#### Frequency ####

function buildFreq!(model,params)

    params.total_freq = zeros(length(params.ebins))
    en_spacing = params.ebins[2]-params.ebins[1]
    for i in 1:length(params.total_freq)
        params.total_freq[i] = maximum([sum([cs(en)*sqrt(2*en*q/me) for (cs,type) in splatCS(model)]) for en in LinRange(params.ebins[i],params.ebins[i]+en_spacing,100)])*1e-20
        
        #TODO: Might be a problem with infinity here
        #TODO: Spactially dependent freq?
    end
end

@xport function splatCS(model)
    zip((model.gas.elastic...,model.gas.excitation...,model.gas.ionisation...,model.gas.attachment...),((:ELASTIC for i in 1:length(model.gas.elastic))...,(:EXCITATION for i in 1:length(model.gas.excitation))...,(:IONIZATION for i in 1:length(model.gas.ionisation))...,(:ATTACHMENT for i in 1:length(model.gas.attachment))...))
end

#### Density ####

function buildDensity!(n0::UniformDensity,params)
    params.n0 = n0.val
end

#### Electric Field ####

function buildField!(n0::UniformDensity,En0::UniformField,params)
    params.a = SVector(0.,0.,q*Td*En0.val*n0.val/me)
end



#### PARTICLE ####

function MCModels.PART(model,params)

    particle = PART()
    particle.r = initPartr(model.r_0,params)
    particle.v = initPartv(model.v_0,params)
    particle.ebin = findMinInd(params.ebins,getEnergy(particle))

    return particle

end

function getEnergy(particle)
    0.5*me*dot(particle.v,particle.v)/q
end

function getVelocity(ε)
    sqrt(ε*2*q/me)
end


function initPartr(r_0::ZeroSpacial,params)
    [0.,0.,0.]
end
function initPartv(v_0::ZeroVelocity,params)
    [0.,0.,0.]
end

function initPartv(v_0::InitVelocity,params)
    v_0.axis .* getVelocity(v_0.ε)
end

end
