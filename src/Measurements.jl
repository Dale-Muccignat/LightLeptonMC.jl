module Measurements

using AndExport, PhysicalConstants.CODATA2018, Unitful, Statistics
using ..Structs

const amu = ustrip(float(AtomicMassConstant))
const q = ustrip(float(ElementaryCharge))
const me = ustrip(float(ElectronMass))
const k = ustrip(float(BoltzmannConstant))
const ħ = ustrip(float(ReducedPlanckConstant))
const Td = 1e-21

@xport function process(type::SWARM,meas,params,model)
    weights!(meas)
    swarmMeas = swarmMeasurements(meas,params,model)
    meas, swarmMeas
end

@xport function process(type::LμJ,meas,params,model)
    meas.measLμJ.EELS ./= model.n
    meas, SWARMMEAS()
end

function weights!(meas)
    for i in 1:length(meas.meast.w)
        meas.meast.εᵢ[i] *= 0.5*me/meas.meast.w[i]/q
        meas.meast.v²[i,:] ./= meas.meast.w[i] 
        meas.meast.r²[i,:] ./= meas.meast.w[i] 
        meas.meast.r[i,:] ./= meas.meast.w[i] 
        meas.meast.v[i,:] ./= meas.meast.w[i]
        meas.meast.rv[i,:] ./= meas.meast.w[i]
    end
end

function swarmMeasurements(meas,params,model)

    bεs = []
    bWs = []
    bDts = []
    bDls = []
    bins = length(meas.meast.w[end÷2:end])
    for i in length(meas.meast.w)-bins:length(meas.meast.w)
        push!(bεs,meas.meast.εᵢ[i])
        push!(bWs,(meas.meast.r[i,3]-meas.meast.r[i-1,3])/params.tbins[2])

		push!(bDls, (1/2 * ( (meas.meast.r²[i,3] - meas.meast.r[i,3]^2) - (meas.meast.r²[i-1,3] - meas.meast.r[i-1,3]^2) )) / (params.tbins[2]))
		push!(bDts, (1/4 * ( (meas.meast.r²[i,2] + meas.meast.r²[i,1]) - (meas.meast.r²[i-1,2] + meas.meast.r²[i-1,1]) )) / (params.tbins[2]))
    end

    SWARMMEAS(mean(bεs),mean(bWs),mean(bDts)*params.n0,mean(bDls)*params.n0)
end




end
