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

    ebins = [0.,(10. .^(-3:0.1:1.5))...]
    tbins = collect(LinRange(0,model.t_f,100))

    params = MCPARAMS(
        ebins=ebins,
        tbins=tbins,
        n0=buildDensity(model.n0),
        a=buildField(model.En0),
        total_freq=buildFreq(model,ebins)
        )

    # buildFreq!(model,params) #TODO
    # buildDensity!(model.n0,params)
    # buildField!(model.n0,model.En0,params)
    return params
end

function initMeas(params::MCPARAMS = MCPARAMS())

    MEAS(params)


end

#### Frequency ####

function buildFreq(model,ebins)

    total_freq = zeros(length(ebins))
    en_spacing = ebins[end]-ebins[end-1]
    for i in 1:length(total_freq)-1
        total_freq[i] = maximum([sum([cs(en)*sqrt(2*en*q/me) for (cs,type) in splatCS(model)]) for en in LinRange(ebins[i]+0.000001,ebins[i+1],1000)])*1e-20
        
        # if total_freq[i] == Inf || isnan(total_freq[i]) total_freq[i] = 1e10 end
        #TODO: Might be a problem with infinity here
        #TODO: Spactially dependent freq?
    end
    total_freq[end] = maximum([sum([cs(en)*sqrt(2*en*q/me) for (cs,type) in splatCS(model)]) for en in LinRange(ebins[end],ebins[end]+en_spacing,1000)])*1e-20
    return total_freq
end

@xport function splatCS(model)
    zip((model.gas.elastic...,model.gas.excitation...,model.gas.ionisation...,model.gas.attachment...),((:ELASTIC for i in 1:length(model.gas.elastic))...,(:EXCITATION for i in 1:length(model.gas.excitation))...,(:IONIZATION for i in 1:length(model.gas.ionisation))...,(:ATTACHMENT for i in 1:length(model.gas.attachment))...))
end

#### Density ####

function buildDensity(n0::UniformDensity)
    n0.val
end

#### Electric Field ####

function buildField(En0::UniformField)
    SVector(0.,0.,q*Td*En0.val/me)
end



#### PARTICLE ####

function MCModels.PART(model,params)

    particle = PART()
    particle.r = initPartr(model.r_0,params)
    particle.v = initPartv(model.v_0,params)
    particle.ebin = findMinInd(params.ebins,getEnergy(particle.v))

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
function initPartv(v_0::ZeroVelocity,params)
    [0.,0.,0.]
end

function initPartv(v_0::InitVelocity,params)
    v_0.axis .* getVelocity(v_0.ε)
end

end
