module Measurements

using AndExport, Statistics, LinearAlgebra, MCSimUnits, Measurements
using ..Structs


@xport function process(type::SWARM,meas,params,model)
    weights!(meas,model)
    swarmMeas = swarmMeasurements(meas,params,model)
    meas, swarmMeas
end

@xport function process(type::LμJ,meas,params,model)
    :ε in model.meas ? meas.measLμJ.ε ./= model.n : nothing
    :εχ in model.meas ? meas.measLμJ.εχ ./= model.n : nothing
    :εχ in model.meas ? meas.measLμJ.εχ_reject ./= model.n : nothing
    :εχ_ybound in model.meas ? meas.measLμJ.εχ ./= model.n : nothing
    :εχ_ybound in model.meas ? meas.measLμJ.εχ_reject ./= model.n : nothing
    :χ in model.meas ? meas.measLμJ.χ ./= model.n : nothing
    :zyt in model.meas ? meas.measLμJ.zyt ./= model.n : nothing
    meas, SWARMMEAS()
end

@xport function process(type::CELL,meas,params,model)
    :ε in model.meas ? meas.measCell.ε ./= model.n : nothing
    :x in model.meas ? meas.measCell.x ./= model.n : nothing
    :y in model.meas ? meas.measCell.y ./= model.n : nothing
    :z in model.meas ? meas.measCell.z ./= model.n : nothing
    :tz in model.meas ? meas.measCell.tz ./= model.n : nothing
    :xy in model.meas ? meas.measCell.xy ./= model.n : nothing
    :xy in model.meas ? meas.measCell.xybins ./= model.n : nothing
    :SP in model.meas ? meas.measg.SP_ε /= model.n : nothing
    :SP_z in model.meas ? meas.measg.SP_εs /= model.n : nothing
    :SP in model.meas ? meas.measg.SP_r /= model.n : nothing
    meas.measg.coschi /= model.n
    weights!(meas,model)
    swarmMeas = swarmMeasurements(meas,params,model)
    meas, swarmMeas
end

function weights!(meas,model)
    for i in 1:length(meas.meast.w) if meas.meast.w[i] !== 0.
        meas.meast.εᵢ[i] /= meas.meast.w[i]
        meas.meast.v²[i,:] ./= meas.meast.w[i]
        meas.meast.r²[i,:] ./= meas.meast.w[i]
        meas.meast.r[i,:] ./= meas.meast.w[i]
        meas.meast.v[i,:] ./= meas.meast.w[i]
        meas.meast.rv[i,:] ./= meas.meast.w[i]
        # meas.measc.r[i,:] ./= meas.meast.w[i]
        # meas.measc.v[i,:] ./= meas.meast.w[i]
        # meas.measc.rv[i,:] ./= meas.meast.w[i]
    end end
    for i in 1:length(meas.measc.weight) if meas.measc.weight[i] !== 0. *s
        meas.measc.r_final[i,:] = meas.measc.r[i,:] ./ meas.measc.weight[i]
        meas.measc.r²_final[i,:] = meas.measc.r²[i,:] ./ meas.measc.weight[i]
        meas.measc.v_final[i,:] = meas.measc.v[i,:] ./ meas.measc.weight[i]
        meas.measc.rv_final[i,:] = meas.measc.rv[i,:] ./ meas.measc.weight[i]
        meas.measc.r²v_final[i,:] = meas.measc.r²v[i,:] ./ meas.measc.weight[i]
        meas.measc.rx²vz_final[i] = meas.measc.rx²vz[i] ./ meas.measc.weight[i]
        meas.measc.rzrxvx_final[i] = meas.measc.rzrxvx[i] ./ meas.measc.weight[i]
    end end
end

function swarmMeasurements(meas,params,model)

    εs::Vector{Quantity} = []
    Rs::Vector{Quantity} = []
    bWs::Vector{Quantity} = []
    bDts::Vector{Quantity} = []
    bDls::Vector{Quantity} = []
    fQts::Vector{Quantity} = []
    fQls::Vector{Quantity} = []
    fεs::Vector{Quantity} = []
    fWs::Vector{Quantity} = []
    fDts::Vector{Quantity} = []
    fDls::Vector{Quantity} = []
    #bins = length(meas.meast.w[1:end])
    for i in 2:length(meas.meast.w)
        push!(εs,meas.meast.εᵢ[i])
        push!(Rs,(log(meas.meast.w[i])-log(meas.meast.w[i-1]))/(params.tbins[2]-params.tbins[1]))

        push!(bWs,(meas.meast.r[i,3]-meas.meast.r[i-1,3])/(params.tbins[2]-params.tbins[1]))

        push!(bDls, (1/2 * ( (meas.meast.r²[i,3] - meas.meast.r[i,3]^2) - (meas.meast.r²[i-1,3] - meas.meast.r[i-1,3]^2) )) / (params.tbins[2]-params.tbins[1]))

        push!(bDts, (1/4 * ( (meas.meast.r²[i,2] + meas.meast.r²[i,1]) - (meas.meast.r²[i-1,2] + meas.meast.r²[i-1,1]) )) / (params.tbins[2]-params.tbins[1]))

        if :f in model.meas
            push!(fWs,sqrt(dot(meas.measc.v_final[i,:],meas.measc.v_final[i,:])))
            push!(fDls, meas.measc.rv_final[i,3] - meas.measc.v_final[i,3]*meas.measc.r_final[i,3])
            push!(fDts, 0.5 * (meas.measc.rv_final[i,2]+meas.measc.rv_final[i,1]))
            push!(fQls, meas.measc.r²v_final[i,3]/2-meas.measc.v_final[i,3]*meas.measc.r²_final[i,3]/2-meas.measc.r_final[i,3]*meas.measc.rv_final[i,3]+meas.measc.r_final[i,3]*meas.measc.r_final[i,3]*meas.measc.v_final[i,3])
            push!(fQts, meas.measc.rx²vz_final[i]/6-meas.measc.rzrxvx_final[i]/3-meas.measc.v_final[i,3]*meas.measc.r²_final[i,1]/6+meas.measc.r_final[i,3]*meas.measc.rv_final[i,1]/3)
        end
    end

    εs = uconvert.(eV,    εs)
    Rs = uconvert.(s^-1,    Rs)
    bWs = uconvert.(m/s,    bWs)
    bDts = uconvert.(m^2/s,    bDts)
    bDls = uconvert.(m^2/s,    bDls)
    fεs = uconvert.(eV,    fεs)
    fWs = uconvert.(m/s,    fWs)
    fDts = uconvert.(m^2/s,    fDts)
    fDls = uconvert.(m^2/s,    fDls)
    fQts = uconvert.(m^3/s,    fQts)
    fQls = uconvert.(m^3/s,    fQls)

    if isempty(fWs)
        fWs = bWs
        fDls = bDls
        fDts = bDts
        fQls = [1m^3/s,1m^3/s]
        fQts = [1m^3/s,1m^3/s]
    end

    swarm=SWARMMEAS(measurement(mean(εs[end÷2:end]),std(εs[end÷2:end])),εs,
        measurement(mean(Rs[end÷2:end]./uconvert(m^-3,params.n0)),std(Rs[end÷2:end]./uconvert(m^-3,params.n0))),Rs./uconvert(m^-3,params.n0),
        measurement(mean(bWs[end÷2:end]),std(bWs[end÷2:end])),bWs,
        measurement(mean(fWs[end÷2:end]),std(fWs[end÷2:end])),fWs,
        measurement(mean(bDts[end÷2:end].*uconvert(m^-3,params.n0)),std(bDts[end÷2:end].*uconvert(m^-3,params.n0))),bDts.*uconvert(m^-3,params.n0),
        measurement(mean(fDts[end÷2:end].*uconvert(m^-3,params.n0)),std(fDts[end÷2:end].*uconvert(m^-3,params.n0))),fDts.*uconvert(m^-3,params.n0),
        measurement(mean(bDls[end÷2:end].*uconvert(m^-3,params.n0)),std(bDls[end÷2:end].*uconvert(m^-3,params.n0))),bDls.*uconvert(m^-3,params.n0),
        measurement(mean(fDls[end÷2:end].*uconvert(m^-3,params.n0)),std(fDls[end÷2:end].*uconvert(m^-3,params.n0))),fDls.*uconvert(m^-3,params.n0),
        measurement(mean(fQls[end÷2:end].*uconvert(m^-3,params.n0).^2),std(fQls[end÷2:end].*uconvert(m^-3,params.n0).^2)),fQls.*uconvert(m^-3,params.n0).^2,
        measurement(mean(fQts[end÷2:end].*uconvert(m^-3,params.n0).^2),std(fQts[end÷2:end].*uconvert(m^-3,params.n0).^2)),fQts.*uconvert(m^-3,params.n0).^2,
    )
    println("""
    Energy: $(swarm.ε)
    Rate: $(swarm.R)
    Bulk:
    W: $(swarm.bW)
    n₀Dₗ: $(swarm.bDl)
    n₀Dₜ: $(swarm.bDt)
    Flux:
    W: $(swarm.fW)
    n₀Dₗ: $(swarm.fDl)
    n₀Dₜ: $(swarm.fDt)
    n₀²Qₗ: $(swarm.fQl)
    n₀²Qₜ: $(swarm.fQt)
        """)
    return swarm
end

end #Module
