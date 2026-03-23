#TODO: Add TimerOutputs
#TODO: Only select electrons with position and velocity directed toward detector in LuJ
#TODO: Tidy up measurement types...

function step_particle!(p, model, params, meas, new_particles; coll_limit = Inf)
    while p.alive && sum(meas.measg.coll_count) < coll_limit
        dt = getDt(p, model, params)

        if :f in model.meas
            measureCumulative!(dt.time, p, model.En0, model.n0, model, params, meas)
        end

        updateParticle!(dt.time, p, model.En0, model.n0, model, params, meas)

        event!(dt, p, model, params, meas, new_particles)

        p.alive &= checkAlive(model.type, p, model, meas, params)

        measureEvent!(model.type, dt, p, model, params, meas)

        checkBins(dt, p, params)
    end

    return measureDeath!(model.type, p, model, params, meas)
end

function flight(particles::Vector{P}, model, params) where {P <: PART}
    :f in model.meas &&
        (params.a[1] > 0m^4 / s^2 || params.a[2] > 0m^4 / s^2) &&
        @warn "Measuring flux with x or y field, third order will be wrong."

    meas = initMeas(model, params)

    for p in particles
        measureTbin!(p, model, params, meas)
    end

    new_particles = copy(particles)

    while !isempty(new_particles)
        p = popfirst!(new_particles)
        step_particle!(p, model, params, meas, new_particles)
    end

    return meas
end

function flightGetParticles(
        particles::Vector{P}, model, params;
        coll_limit = Inf, t_limit = Inf * s
    ) where {P <: PART}

    :f in model.meas &&
        (params.a[1] > 0m^4 / s^2 || params.a[2] > 0m^4 / s^2) &&
        @warn "Measuring flux with x or y field, third order will be wrong."

    meas = initMeas(model, params)

    measureTbin!(particles, model, params, meas)

    new_particles = particles
    i = 1

    while i ≤ length(new_particles)
        if length(new_particles) - i > 50000
            error("Too many ionisation events")
        end

        p = new_particles[i]

        step_particle!(
            p, model, params, meas, new_particles;
            coll_limit = coll_limit
        )

        i += 1
    end

    return filter!(p -> p.weight != 0.0, new_particles)
end

function checkBins(dt, particle, params)
    # ebin is not a tracked quantity
    # if params.ebins |> !isempty && particle.ebin < length(params.ebins)
    #     (params.ebins[particle.ebin]-(eps()+1e-6)*eV) <= getEnergy(particle.v) <= (params.ebins[particle.ebin+1]+(eps()+1e-6)*eV) || error("""ebin error
    #         $(uconvert(eV,params.ebins[particle.ebin]))
    #         $(uconvert(eV,getEnergy(particle.v)))
    #         $(uconvert(eV,params.ebins[particle.ebin+1]))
    #         $(dt)
    #                                                                                                                                """)
    # end
    if params.cbins |> !isempty && particle.cbin < length(params.cbins)
        (params.cbins[particle.cbin] - (eps() + 1.0e-8) * eV) <= getEnergy(particle.v) <= (params.cbins[particle.cbin + 1] + (eps() + 1.0e-8) * eV) || error(
            """CBin error
            $(uconvert(eV, params.cbins[particle.cbin]))
            $(uconvert(eV, getEnergy(particle.v)))
            $(uconvert(eV, params.cbins[particle.cbin + 1]))
            $(dt)
            """
        )
    end
    if params.zbins |> !isempty && 0 < particle.zbin < length(params.zbins)
        if (params.zbins[particle.zbin] - (1.0e-12 + eps()) * m̃) <= particle.r[3] <= (params.zbins[particle.zbin + 1] + (1.0e-12 + eps()) * m̃) && particle.r[3] >= params.zbins[1]
        else
            error(
                """
                zbin error
                $(params.zbins[1])
                $(params.zbins[particle.zbin] - (1.0e-12 + eps()) * m̃)
                $(params.zbins[particle.zbin])
                $(particle.r[3])
                $(params.zbins[particle.zbin + 1])
                $(params.zbins[particle.zbin + 1] + (1.0e-12 + eps()) * m̃)
                $(dt)
                """
            )
        end
    end
    return if params.tbins |> !isempty && particle.tbin < length(params.tbins) && particle.t > params.tbins[1]
        (params.tbins[particle.tbin] - eps() * s) <= particle.t <= (params.tbins[particle.tbin + 1] + eps() * s) || error(
            """tbin error
            $(params.tbins[particle.tbin])
            $(particle.t)
            $(params.tbins[particle.tbin + 1])
            $(dt)
            """
        )
        #= $(debugParticle("",particle_states[1])) =#
        #= $(debugDts(dts_states[1])) =#
        #= $(debugParticle("",particle_states[2])) =#
        #= $(debugDts(dts_states[2])) =#
        #= $(debugParticle("",particle_states[3])) =#
        #= $(debugDts(dts_states[3])) =#
    end
end

function checkAlive(type, particle, model, meas_inner, params)
    alive = true
    if params.tbins |> !isempty
        alive = alive && particle.t < model.type.t_f #*params.time_scale
        alive || @debug "Particle timed out"
        # alive = alive && sum(meas_inner.measg.coll_count) < 5
    end

    return alive && checkAliveModel(type, particle, model, meas_inner, params)
end

function checkAliveModel(type::SWARM, particle, model, meas_inner, params)
    return true
end

function checkAliveModel(type::LμJ, particle, model, meas_inner, params)
    if getEnergy(particle.v) < type.thresh #*params.energy_scale
        @debug "Absorbed"
    end

    if particle.rbin == length(params.rbins)
        @debug "Detected"
    end

    return particle.rbin != length(params.rbins) && getEnergy(particle.v) > type.thresh #&& sumColl(meas_inner.measg.coll_count) < 2
end

function checkAliveModel(type::CELL, particle, model, meas_inner, params)
    alive = true
    # if particle.r[3] > 15mm
    #     alive = false
    # end

    if getEnergy(particle.v) < type.thresh #*params.energy_scale
        @debug "Absorbed"
    end

    if params.zbins |> !isempty
        # :SP_z in model.meas ? alive = alive && particle.zbin != length(params.zbins) : alive = alive && particle.zbin != length(params.zbins) && particle.zbin != 0
        #= alive = alive && particle.zbin != length(params.zbins) && particle.zbin > 1 =#
        alive = alive && particle.zbin != length(params.zbins) && particle.zbin != 1
        (particle.zbin == length(params.zbins) || particle.zbin == 1) && @debug "Detected Z"
    end

    if params.ybins |> !isempty
        alive = alive && particle.ybin != length(params.ybins)
        particle.ybin == length(params.ybins) && @debug "Detected Y"
    end

    if params.xbins |> !isempty
        alive = alive && particle.xbin != length(params.xbins)
        particle.xbin == length(params.xbins) && @debug "Detected X"
    end

    :SP in model.meas ? alive = alive && sum(meas_inner.measg.coll_count) == 0 : nothing

    return alive && getEnergy(particle.v) > type.thresh #*params.energy_scale
end

function measureEvent!(type::CELL, dt, particle, model, params, meas_inner)
    if dt.type == :dt_zbin && dt.direction == :Up
        if :z in model.meas
            if meas_inner.measCell.z[particle.zbin] == 0.0
                meas_inner.measCell.z[particle.zbin] += particle.weight
            end
        end
        if :εz in model.meas
            # if sum(meas_inner.measCell.εz[:,particle.zbin]) == 0
            # meas_inner.measCell.εz[particle.ebin,particle.zbin] += uconvert(eV,getEnergy(particle.v))
            # meas_inner.measCell.εz[particle.ebin,particle.zbin] += particle.weight
            # end
        end
    end
    return if dt.type == :dt_zbin
        if :εz in model.meas && particle.zbin != 0
            meas_inner.measCell.εz_w[particle.zbin] += particle.weight
            meas_inner.measCell.εz[particle.zbin] += particle.weight * uconvert(eV, getEnergy(particle.v))
        end
    end
end

function measureEvent!(type::LμJ, dt, particle, model, params, meas_inner)
end
function measureEvent!(type::SWARM, dt, particle, model, params, meas_inner)
end

@fastmath function measureCumulative!(dt, particle, En0, n0, model, params, meas_inner)
    # e_bin += particle.weight * 0.5*consts.me*(dot(particle.v,particle.v)+2*sum(particle.v.*vars.a)*dt+dot(vars.a,vars.a)*dt^2)*dt
    a = params.a * params.n0
    r = particle.r
    v = particle.v
    tbin = particle.tbin
    w = particle.weight
    @inbounds @views for i in 1:3
        meas_inner.measc.r[tbin, i] += w * (r[i] + 0.5 * v[i] * dt + (1 / 6) * a[i] * (dt^2)) * dt
        meas_inner.measc.r²[tbin, i] += w * dt * (
            a[i]^2 * dt^4 / 20 +
                dt^2 * (a[i] .* r[i] + v[i] .^ 2) / 3 +
                a[i] * dt^3 .* v[i] / 4 +
                r[i]^2 +
                r[i] * v[i] * dt
        )

        meas_inner.measc.v[tbin, i] += w * (v[i] + 0.5 * a[i] * dt) * dt
        meas_inner.measc.rv[tbin, i] += w * (r[i] * v[i] + 0.5 * (r[i] * a[i] + v[i]^2) * dt + 0.5 * v[i] * a[i] * (dt^2) + (1 / 8) * (a[i]^2) * (dt^3)) * dt
        # NOTE: These have been optimised for electric field in z direction
        #
        meas_inner.measc.r²v[tbin, i] += 1 / 24 * dt * w * (
            a[i] * dt + 2 * v[i]
        ) * (
            6 * r[i] * dt * (a[i] * dt + 2 * v[i]) + dt^2 * (a[i] * dt + 2v[i])^2 + 12r[i]^2
        )

    end
    @inbounds begin
        meas_inner.measc.rx²vz[tbin] += 1 / 12 * dt^2 * w * (
            4 * dt * v[1] * (2 * a[3] * r[1] + v[3] * v[1]) +
                6 * r[1] * (a[3] * r[1] + 2 * v[1] * v[3]) +
                3 * a[3] * dt^2 * v[1]^2 +
                12 * v[3] * r[1]^2 / dt
        )

        meas_inner.measc.rzrxvx[tbin] += w * v[1] * dt * (
            a[3] * v[1] * dt^3 / 8 +
                dt^2 * (a[3] * r[1] + 2 * v[1] * v[3]) / 6 +
                dt * (v[1] * r[3] + v[3] * r[1]) / 2 +
                r[3] * r[1]
        )

        meas_inner.measc.weight[tbin] += w * dt
        meas_inner.measc.dt[tbin] += dt
    end
end

function measureDeath!(type::SWARM, particle, model, params, meas_inner)
end

function measureDeath!(type::LμJ, particle, model, params, meas_inner)
    return if particle.rbin == length(params.rbins) && sum(meas_inner.measg.coll_count) != 0
        :χ in model.meas ? meas_inner.measLμJ.χ[findMinInd(params.χbins, abs(atan(rad, particle.r[2], particle.r[3])))] += particle.weight : nothing

        if :εχ in model.meas
            @debug """
            Velocity angle:  $(abs(atan(°, particle.v[2], particle.v[3])))
            Spatial angle:  $(abs(atan(°, particle.r[2], particle.r[3])))
            Velocity angle:  $(atan(°, particle.v[2], particle.v[3]))
            Spatial angle:  $(atan(°, particle.r[2], particle.r[3]))
            """
            angle_position = abs(atan(°, particle.r[2], particle.r[3]))
            angle_velocity = abs(atan(°, particle.v[2], particle.v[3]))
            if abs(angle_position - angle_velocity) < type.Δχ_detector #&& sumColl(meas_inner.measg.coll_count) == 1
                meas_inner.measLμJ.εχ[findMinInd(params.ebins, getEnergy(particle.v)), findMinInd(params.χbins, angle_position)] += particle.weight
                @debug "accepted"
            else
                meas_inner.measLμJ.εχ_reject[findMinInd(params.ebins, getEnergy(particle.v)), findMinInd(params.χbins, angle_position)] += particle.weight
                @debug "rejected"
            end
        elseif :ε in model.meas
            angle_position = abs(atan(°, particle.r[2], particle.r[3]))
            angle_velocity = abs(atan(°, particle.v[2], particle.v[3]))
            if abs(angle_position - angle_velocity) < type.Δχ_detector #&& sumColl(meas_inner.measg.coll_count) == 1
                meas_inner.measLμJ.ε[findMinInd(params.ebins, getEnergy(particle.v))] += particle.weight
            end
        elseif :εχ_ in model.meas
            angle_position = abs(atan(°, particle.r[2], particle.r[3]))
            angle_velocity = abs(atan(°, particle.v[2], particle.v[3]))
            height_position = abs(particle.r[1])
            @debug """
            Velocity angle:  $(abs(atan(°, particle.v[2], particle.v[3])))
            Spatial angle:  $(abs(atan(°, particle.r[2], particle.r[3])))
            Velocity angle:  $(atan(°, particle.v[2], particle.v[3]))
            Spatial angle:  $(atan(°, particle.r[2], particle.r[3]))
            height: $(height_position)
            """
            if abs(angle_position - angle_velocity) < type.Δχ_detector #&& height_position < 5mm
                # TODO: Add detection at +- 1eV
                inds = findMinInd(params.ebins, getEnergy(particle.v) - type.Δε_detector):findMinInd(params.ebins, getEnergy(particle.v) + type.Δε_detector)
                inds_angle = findMinInd(params.χbins, angle_position - type.Δχ_detector / 2):findMinInd(params.χbins, angle_position + type.Δχ_detector / 2)
                meas_inner.measLμJ.εχ[inds, inds_angle] .+= particle.weight
                @debug "accepted - $(inds) $(inds_angle)"
            else
                inds = findMinInd(params.ebins, getEnergy(particle.v) - type.Δε_detector / 2):findMinInd(params.ebins, getEnergy(particle.v) + type.Δε_detector / 2)
                inds_angle = findMinInd(params.χbins, angle_position - type.Δχ_detector / 2):findMinInd(params.χbins, angle_position + type.Δχ_detector / 2)
                meas_inner.measLμJ.εχ[inds, inds_angle] .+= particle.weight
                @debug "rejected"
            end
        end


    end
end

function measureDeath!(type::CELL, particle, model, params, meas_inner)
    if :z in model.meas && particle.zbin == length(params.zbins)
        :xy in model.meas ? meas_inner.measCell.xy[findMinInd(meas_inner.measCell.xybins, sqrt(particle.r[1]^2 + particle.r[2]^2))] += particle.weight : nothing
        :ε in model.meas ? meas_inner.measCell.ε[findMinInd(params.ebins, getEnergy(particle.v))] += particle.weight : nothing
    end
    return if :tz in model.meas && particle.zbin == length(params.zbins)
        # @show "here"
        # @show particle.weight*(particle.t-params.tbins[particle.tbin])/(params.tbins[2]-params.tbins[1])
        # @show uconvert(s,particle.t), particle.weight, particle.zbin, particle.tbin
        # @show meas_inner.measCell.tz[particle.tbin,particle.zbin-1]
        # meas_inner.measCell.tz[particle.tbin,particle.zbin-1] -= particle.weight*(particle.t-params.tbins[particle.tbin])/(params.tbins[2]-params.tbins[1])
        # meas_inner.measCell.tz[particle.tbin+1,particle.zbin-1] += particle.weight*(particle.t-params.tbins[particle.tbin])/(params.tbins[2]-params.tbins[1])
        # @show meas_inner.measCell.tz[particle.tbin,particle.zbin-1]
        # meas_inner.measCell.v_tz[particle.tbin,particle.zbin-1] += particle.weight*particle.v[3]*(particle.t-params.tbins[particle.tbin])/(params.tbins[2]-params.tbins[1])
    end
end

function measureCollision!(energy_loss, particle, model, params, meas_inner)
    return if :SP in model.meas && sum(meas_inner.measg.coll_count) == 1
        meas_inner.measg.SP_ε = energy_loss
        meas_inner.measg.SP_r = sqrt(dot(particle.r, particle.r)) # Assumes particle started at [0,0,0]
    end
end

function event!(dt, particle, model, params, meas_inner, new_particles)
    return if dt.type == :dt_collision
        eventCollision!(dt, particle, model, params, meas_inner, new_particles)
    elseif dt.type == :dt_energy
        eventEnergy!(dt, particle, model, params, meas_inner)
    elseif dt.type == :dt_tbin
        eventTbin!(dt, particle, model, params, meas_inner)
    elseif dt.type == :dt_rbin
        eventRBin!(dt, particle, model, params, meas_inner)
    elseif dt.type == :dt_zbin
        eventZBin!(dt, particle, model, params, meas_inner)
    elseif dt.type == :dt_ybin
        eventYBin!(dt, particle, model, params, meas_inner)
    else
    end
end

@fastmath function updateParticle!(dt, particle, En0::UniformField, n0::UniformDensity, model, params, meas_inner)
    particle.r += (particle.v + 0.5 * params.a * params.n0 * dt) * dt
    particle.v += params.a * params.n0 * dt
    particle.t += dt
end

@fastmath function updateParticle!(dt, particle, En0::ZeroField, n0, model, params, meas_inner)
    particle.r += particle.v * dt
    particle.t += dt
end
#### * TBins Measurements
function eventTbin!(dt, particle, model, params, meas_inner)
    #' @debug "TBin Up: #$(particle.zbin) $(dt.time)"
    particle.tbin += 1
    return measureTbin!(particle, model, params, meas_inner)
end

@inbounds function measureTbin!(particle, model, params, meas_inner)
    r = particle.r
    v = particle.v
    w = particle.weight
    if :t ∈ model.meas
        @views for i in 1:3
            meas_inner.meast.v²[particle.tbin, i] += v[i]^2 * w
            meas_inner.meast.r²[particle.tbin, i] += r[i]^2 * w
            meas_inner.meast.r³[particle.tbin, i] += r[i]^3 * w
            meas_inner.meast.r[particle.tbin, i] += r[i] * w
            meas_inner.meast.v[particle.tbin, i] += v[i] * w
            meas_inner.meast.rv[particle.tbin, i] += v[i] * r[i] * w
        end
        meas_inner.meast.zx²[particle.tbin] += particle.r[3] * particle.r[2]^2 * w
        meas_inner.meast.w[particle.tbin] += w
        meas_inner.meast.εᵢ[particle.tbin] += 0.5 * me * dot(v, v) * w
    end

    if :zyt ∈ model.meas
        meas_inner.measLμJ.zyt[particle.zbin, particle.ybin, particle.tbin] += w
    end

    if :tz ∈ model.meas
        meas_inner.measCell.tz[particle.tbin, particle.zbin] += w
    end
    if :tz_v ∈ model.meas
        meas_inner.measCell.v_tz[particle.tbin, particle.zbin] += w * v[3]
    end
    # :tz ∈ model.meas && params.zbins[2] < particle.r[3] < params.zbins[end] ? meas_inner.measCell.tz[particle.tbin,1] += w : nothing
    # :tz ∈ model.meas && particle.r[3] < params.zbins[end] ? meas_inner.measCell.tz[particle.tbin,1] += w : nothing
end


#### * Update Bin numbers

function updateCbin!(particle, params)
    ea = getEnergy(particle.v)
    while ea < params.cbins[particle.cbin]
        particle.cbin -= 1
    end
    return if particle.cbin != length(params.cbins)
        while ea > params.cbins[particle.cbin + 1]
            particle.cbin += 1
        end
    end
    # if params.ebins !== nothing
    #     while ea < params.ebins[particle.ebin] && particle.ebin != 1
    #         particle.ebin-= 1
    #     end
    #     # if particle.ebin != length(params.ebins)
    #         while particle.ebin != length(params.ebins) && ea > params.ebins[particle.ebin+1]
    #             particle.ebin+= 1
    #         end
    #     # end
    # end
end
# * Fetch Functions
