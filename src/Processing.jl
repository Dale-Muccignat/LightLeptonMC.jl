function process(type::SWARM, meas, params, model)
    meas = weights(meas, model)
    swarmMeas = swarmMeasurements(meas, params, model)
    meas, swarmMeas
end

function process(type::LμJ, meas, params, model)
    #= :ε in model.meas ? meas.measLμJ.ε ./= model.n : nothing =#
    #= :εχ in model.meas ? meas.measLμJ.εχ ./= model.n : nothing =#
    #= :εχ in model.meas ? meas.measLμJ.εχ_reject ./= model.n : nothing =#
    #= :εχ_ybound in model.meas ? meas.measLμJ.εχ ./= model.n : nothing =#
    #= :εχ_ybound in model.meas ? meas.measLμJ.εχ_reject ./= model.n : nothing =#
    #= :χ in model.meas ? meas.measLμJ.χ ./= model.n : nothing =#
    #= :zyt in model.meas ? meas.measLμJ.zyt ./= model.n : nothing =#
    meas, SWARMMEAS()
end

function process(type::CELL, meas, params, model)
    #= :ε in model.meas ? meas.measCell.ε ./= model.n : nothing =#
    #= :x in model.meas ? meas.measCell.x ./= model.n : nothing =#
    #= :y in model.meas ? meas.measCell.y ./= model.n : nothing =#
    #= :z in model.meas ? meas.measCell.z ./= model.n : nothing =#
    #= :tz in model.meas ? meas.measCell.tz ./= model.n : nothing =#
    #= :xy in model.meas ? meas.measCell.xy ./= model.n : nothing =#
    #= :xy in model.meas ? meas.measCell.xybins ./= model.n : nothing =#
    #= :SP in model.meas ? meas.measg.SP_ε /= model.n : nothing =#
    #= :SP_z in model.meas ? meas.measg.SP_εs /= model.n : nothing =#
    #= :SP in model.meas ? meas.measg.SP_r /= model.n : nothing =#
    #= meas.measg.coschi /= model.n =#

    # weights!(meas,model)
    # swarmMeas = swarmMeasurements(meas,params,model)
    # meas, swarmMeas
    meas, SWARMMEAS()
end

function weights(meas, model)
    for i in 1:length(meas.meast.w)
        if meas.meast.w[i] !== 0.0
            meas.meast.εᵢ[i] /= meas.meast.w[i]
            meas.meast.v²[i, :] ./= meas.meast.w[i]
            meas.meast.r²[i, :] ./= meas.meast.w[i]
            meas.meast.r³[i, :] ./= meas.meast.w[i]
            meas.meast.r[i, :] ./= meas.meast.w[i]
            meas.meast.v[i, :] ./= meas.meast.w[i]
            meas.meast.rv[i, :] ./= meas.meast.w[i]
            # meas.measc.r[i,:] ./= meas.meast.w[i]
            # meas.measc.v[i,:] ./= meas.meast.w[i]
            # meas.measc.rv[i,:] ./= meas.meast.w[i]
        end
    end
    for i in 1:length(meas.measc.weight)
        if meas.measc.weight[i] !== 0.0 * s
            meas.measc.r_final[i, :] = meas.measc.r[i, :] ./ meas.measc.weight[i]
            meas.measc.r²_final[i, :] = meas.measc.r²[i, :] ./ meas.measc.weight[i]
            meas.measc.v_final[i, :] = meas.measc.v[i, :] ./ meas.measc.weight[i]
            meas.measc.rv_final[i, :] = meas.measc.rv[i, :] ./ meas.measc.weight[i]
            meas.measc.r²v_final[i, :] = meas.measc.r²v[i, :] ./ meas.measc.weight[i]
            meas.measc.rx²vz_final[i] = meas.measc.rx²vz[i] ./ meas.measc.weight[i]
            meas.measc.rzrxvx_final[i] = meas.measc.rzrxvx[i] ./ meas.measc.weight[i]
        end
    end
    return meas
end

function unweights(meas_scaled, model)
    for i in 1:length(meas_scaled.meast.w)
        if meas_scaled.meast.w[i] !== 0.0
            meas_scaled.meast.εᵢ[i] *= meas_scaled.meast.w[i]
            meas_scaled.meast.v²[i, :] .*= meas_scaled.meast.w[i]
            meas_scaled.meast.r²[i, :] .*= meas_scaled.meast.w[i]
            meas_scaled.meast.r³[i, :] .*= meas_scaled.meast.w[i]
            meas_scaled.meast.r[i, :] .*= meas_scaled.meast.w[i]
            meas_scaled.meast.v[i, :] .*= meas_scaled.meast.w[i]
            meas_scaled.meast.rv[i, :] .*= meas_scaled.meast.w[i]
            # meas_scaled.measc.r[i,:] .*= meas_scaled.meast.w[i]
            # meas_scaled.measc.v[i,:] .*= meas_scaled.meast.w[i]
            # meas_scaled.measc.rv[i,:] .*= meas_scaled.meast.w[i]
        end
    end
    return meas_scaled
end

function swarmMeasurements(meas, params, model)

    εs = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    Rs = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    bWs = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    bDts = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    bDls = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    bQts = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    bQls = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    fQts = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    fQls = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    fεs = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    fWs = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    fDts = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    fDls = Array{Quantity,1}(undef, length(meas.meast.w) - 1)
    #bins = length(meas.meast.w[1:end])
    for i in 2:length(meas.meast.w)
        εs[i-1] = meas.meast.εᵢ[i]
        Rs[i-1] = (log(meas.meast.w[i]) - log(meas.meast.w[i-1])) / (params.tbins[2] - params.tbins[1])

        bWs[i-1] = (meas.meast.r[i, 3] -
                    meas.meast.r[i-1, 3]) /
                   (params.tbins[2] - params.tbins[1])

        bDls[i-1] = (1 / 2 * (
            (meas.meast.r²[i, 3] - meas.meast.r[i, 3]^2) -
            (meas.meast.r²[i-1, 3] - meas.meast.r[i-1, 3]^2))) /
                    (params.tbins[2] - params.tbins[1])

        bDts[i-1] = (1 / 4 * (
            (meas.meast.r²[i, 2] + meas.meast.r²[i, 1]) -
            (meas.meast.r²[i-1, 2] + meas.meast.r²[i-1, 1]))) /
                    (params.tbins[2] - params.tbins[1])

        bQls[i-1] = (1 / 6 * (
            (meas.meast.r³[i, 3] - 3 * meas.meast.r[i, 3] * meas.meast.r²[i, 3] + 2 * meas.meast.r[i, 3]^3) -
            (meas.meast.r³[i-1, 3] - 3 * meas.meast.r[i-1, 3] * meas.meast.r²[i-1, 3] + 2 * meas.meast.r[i-1, 3]^3))) /
                    (params.tbins[2] - params.tbins[1])

        bQts[i-1] = (1 / 6 * (
            (meas.meast.zx²[i] - meas.meast.r[i, 3] * meas.meast.r²[i, 2]) -
            (meas.meast.zx²[i-1] - meas.meast.r[i-1, 3] * meas.meast.r²[i-1, 2]))) /
                    (params.tbins[2] - params.tbins[1])

        if :f in model.meas
            fWs[i-1] = sqrt(dot(meas.measc.v_final[i, :], meas.measc.v_final[i, :]))
            fDls[i-1] = meas.measc.rv_final[i, 3] - meas.measc.v_final[i, 3] * meas.measc.r_final[i, 3]
            fDts[i-1] = 0.5 * (meas.measc.rv_final[i, 2] + meas.measc.rv_final[i, 1])
            fQls[i-1] = meas.measc.r²v_final[i, 3] / 2 - meas.measc.v_final[i, 3] * meas.measc.r²_final[i, 3] / 2 - meas.measc.r_final[i, 3] * meas.measc.rv_final[i, 3] + meas.measc.r_final[i, 3] * meas.measc.r_final[i, 3] * meas.measc.v_final[i, 3]
            fQts[i-1] = meas.measc.rx²vz_final[i] / 6 - meas.measc.rzrxvx_final[i] / 3 - meas.measc.v_final[i, 3] * meas.measc.r²_final[i, 1] / 6 + meas.measc.r_final[i, 3] * meas.measc.rv_final[i, 1] / 3
        end
    end


    if !isassigned(fWs, 1)
        fWs = bWs
        fDls = bDls
        fDts = bDts
        fQls = [1m^3 / s, 1m^3 / s]
        fQts = [1m^3 / s, 1m^3 / s]
    end

    εs = uconvert.(eV, εs)
    Rs = uconvert.(s^-1, Rs)
    bWs = uconvert.(m / s, bWs)
    bDts = uconvert.(m^2 / s, bDts)
    bDls = uconvert.(m^2 / s, bDls)
    # fεs = uconvert.(eV,    fεs)
    fWs = uconvert.(m / s, fWs)
    fDts = uconvert.(m^2 / s, fDts)
    fDls = uconvert.(m^2 / s, fDls)
    bQts = uconvert.(m^3 / s, bQts)
    bQls = uconvert.(m^3 / s, bQls)
    fQts = uconvert.(m^3 / s, fQts)
    fQls = uconvert.(m^3 / s, fQls)

    swarm = SWARMMEAS(measurement(mean(εs[end÷2:end-1]),
            std(εs[end÷2:end-1])), εs,
        measurement(mean(Rs[end÷2:end-1] ./ uconvert(m^-3, params.n0)),
            std(Rs[end÷2:end-1] ./ uconvert(m^-3, params.n0))), Rs ./ uconvert(m^-3, params.n0),
        measurement(mean(bWs[end÷2:end-1]),
            std(bWs[end÷2:end-1])), bWs,
        measurement(mean(fWs[end÷2:end-1]),
            std(fWs[end÷2:end-1])), fWs,
        measurement(mean(bDts[end÷2:end-1] .* uconvert(m^-3, params.n0)),
            std(bDts[end÷2:end-1] .* uconvert(m^-3, params.n0))), bDts .* uconvert(m^-3, params.n0),
        measurement(mean(fDts[end÷2:end-1] .* uconvert(m^-3, params.n0)),
            std(fDts[end÷2:end-1] .* uconvert(m^-3, params.n0))), fDts .* uconvert(m^-3, params.n0),
        measurement(mean(bDls[end÷2:end-1] .* uconvert(m^-3, params.n0)),
            std(bDls[end÷2:end-1] .* uconvert(m^-3, params.n0))), bDls .* uconvert(m^-3, params.n0),
        measurement(mean(fDls[end÷2:end-1] .* uconvert(m^-3, params.n0)),
            std(fDls[end÷2:end-1] .* uconvert(m^-3, params.n0))), fDls .* uconvert(m^-3, params.n0),
        measurement(mean(bQls[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2),
            std(bQls[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2)), bQls .* uconvert(m^-3, params.n0) .^ 2,
        measurement(mean(bQts[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2),
            std(bQts[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2)), bQts .* uconvert(m^-3, params.n0) .^ 2,
        measurement(mean(fQls[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2),
            std(fQls[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2)), fQls .* uconvert(m^-3, params.n0) .^ 2,
        measurement(mean(fQts[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2),
            std(fQts[end÷2:end-1] .* uconvert(m^-3, params.n0) .^ 2)), fQts .* uconvert(m^-3, params.n0) .^ 2,
    )
    println("""
    Energy: $(swarm.ε)
    Rate: $(swarm.R)
    Bulk:
    W: $(swarm.bW)
    n₀Dₗ: $(swarm.bDl)
    n₀Dₜ: $(swarm.bDt)
    n₀²Qₗ: $(swarm.bQl)
    n₀²Qₜ: $(swarm.bQt)
    Flux:
    W: $(swarm.fW)
    n₀Dₗ: $(swarm.fDl)
    n₀Dₜ: $(swarm.fDt)
    n₀²Qₗ: $(swarm.fQl)
    n₀²Qₜ: $(swarm.fQt)
        """)
    return swarm
end

