module Measurements

using AndExport, PhysicalConstants.CODATA2018, Unitful, Statistics, LinearAlgebra
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
    :ε in model.meas ? meas.measLμJ.ε ./= model.n : nothing
    :εχ in model.meas ? meas.measLμJ.εχ ./= model.n : nothing
    :χ in model.meas ? meas.measLμJ.χ ./= model.n : nothing
    :zyt in model.meas ? meas.measLμJ.zyt ./= model.n : nothing
    meas, SWARMMEAS()
end

@xport function process(type::CELL,meas,params,model)
    :ε in model.meas ? meas.measCell.ε ./= model.n : nothing
    :x in model.meas ? meas.measCell.x ./= model.n : nothing
    :y in model.meas ? meas.measCell.y ./= model.n : nothing
    :z in model.meas ? meas.measCell.z ./= model.n : nothing
    :xy in model.meas ? meas.measCell.xy ./= model.n : nothing
    :SP in model.meas ? meas.measg.SP_ε /= model.n : nothing
    :SP in model.meas ? meas.measg.SP_r /= model.n : nothing
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
        meas.measc.r[i,:] ./= meas.meast.w[i] 
        meas.measc.v[i,:] ./= meas.meast.w[i]
        meas.measc.rv[i,:] ./= meas.meast.w[i]
    end
end

function swarmMeasurements(meas,params,model)

    εs::Vector{Float64} = []
    bWs::Vector{Float64} = []
    bDts::Vector{Float64} = []
    bDls::Vector{Float64} = []
    fεs::Vector{Float64} = []
    fWs::Vector{Float64} = []
    fDts::Vector{Float64} = []
    fDls::Vector{Float64} = []
    bins = length(meas.meast.w[end÷2:end])
    for i in length(meas.meast.w)-bins:length(meas.meast.w)
        push!(εs,meas.meast.εᵢ[i])

        push!(bWs,(meas.meast.r[i,3]-meas.meast.r[i-1,3])/params.tbins[2])

		push!(bDls, (1/2 * ( (meas.meast.r²[i,3] - meas.meast.r[i,3]^2) - (meas.meast.r²[i-1,3] - meas.meast.r[i-1,3]^2) )) / (params.tbins[2]))

		push!(bDts, (1/4 * ( (meas.meast.r²[i,2] + meas.meast.r²[i,1]) - (meas.meast.r²[i-1,2] + meas.meast.r²[i-1,1]) )) / (params.tbins[2]))
        if :f in model.meas
            push!(fWs,sqrt(dot(meas.measc.v[i,:],meas.measc.v[i,:])))
            push!(fDls, meas.measc.rv[i,3] - meas.measc.v[i,3]*meas.measc.r[i,3])
            push!(fDts, 0.5 * (meas.measc.rv[i,2]+meas.measc.rv[i,1]))
        end
        
    end

    SWARMMEAS(mean(εs),
              εs,mean(bWs),
              bWs,mean(fWs),
              fWs,mean(bDts),
              bDts.*params.n0,mean(fDts)*params.n0,
              fDts.*params.n0,mean(bDls)*params.n0,
              bDls.*params.n0,mean(fDls)*params.n0,
              fDls.*params.n0
              )
end




end
