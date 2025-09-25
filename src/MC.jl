#TODO: Add TimerOutputs
#TODO: Only select electrons with position and velocity directed toward detector in LuJ
#TODO: Tidy up measurement types...

function flight(particle_og::PART, model, params)
    #' @debug "Conducting flight..."
    :f in model.meas && (params.a[1] > 0m^4 / s^2 || params.a[2] > 0m^4 / s^2) && @warn "Measuring flux with x or y field, third order will be wrong."
    meas_inner = initMeas(model, params)
    # Preallocate arrays
    #= dts = Array{dt_EVENT,1}(undef, 7) =#
    #= dt_times = Array{Quantity,1}(undef, 7) =#
    dts = repeat([dt_EVENT(1.0 * s̃, :before, :before)], 7)
    dt_times = repeat([1.0 * s̃], 7)

    measureTbin!(particle_og, model, params, meas_inner) ## measure initial stuff

    new_particles = [particle_og]
    while !isempty(new_particles)
        #= @show "start" =#
        particle_inner = new_particles[1]
        length(new_particles) % 1000 == 0 && @info "Particles num: " length(new_particles)
        if length(new_particles) > 50000
            #= error("Too many ionisation events") =#
        end
        event_number = 0
        #= particle_states = fill(deepcopy(particle_inner), 3) =#
        #= dts_states = fill(fill(dt_EVENT(1.0 * s̃, :before, :before), 7), 3) =#
        # @show length(new_particles)
        while particle_inner.alive
            event_number += 1
            #@debug "Particle before: " particle
            #' @debug debugParticle("before", particle_inner)
            #= checkBins(dt_EVENT(1.0 * s̃, :before, :before), particle_inner, params, particle_states, dts_states) =#
            checkBins(dt_EVENT(1.0 * s̃, :before, :before), particle_inner, params)
            dt = getDt!(dts, dt_times, particle_inner, model, params)

            :f in model.meas ? measureCumulative!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner) : nothing
            updateParticle!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner)
            #' @debug debugParticle("after flight", particle_inner)

            event!(dt, particle_inner, model, params, meas_inner, new_particles)
            #= particle_states[event_number%3+1] = deepcopy(particle_inner) =#
            #= dts_states[event_number%3+1] = dts =#

            #' @debug debugParticle("after event", particle_inner)
            #= @debug sum(meas_inner.measg.coll_count) > 5 && break =#

            particle_inner.alive = particle_inner.alive && checkAlive(model.type, particle_inner, model, meas_inner, params)
            measureEvent!(model.type, dt, particle_inner, model, params, meas_inner)
            # if verbose && meas_inner.measg.coll_count.elastic > 5 alive = false end
            #= checkBins(dt, particle_inner, params, particle_states, dts_states) =#
            checkBins(dt, particle_inner, params)
        end
        measureDeath!(model.type, particle_inner, model, params, meas_inner)
        popfirst!(new_particles)
    end
    return meas_inner
end

function flight(particle_og::Vector{<:PART}, model, params)
    #' @debug "Conducting flight..."
    :f in model.meas && (params.a[1] > 0m^4 / s^2 || params.a[2] > 0m^4 / s^2) && @warn "Measuring flux with x or y field, third order will be wrong."
    meas_inner = initMeas(model, params)
    # Preallocate arrays
    #= dts = Array{dt_EVENT,1}(undef, 7) =#
    #= dt_times = Array{Quantity,1}(undef, 7) =#
    dts = repeat([dt_EVENT(1.0 * s̃, :before, :before)], 7)
    dt_times = repeat([1.0 * s̃], 7)
    og_num = length(particle_og)

    for p in particle_og
        measureTbin!(p, model, params, meas_inner) ## measure initial stuff
    end

    new_particles = particle_og
    particle_no = 0
    while !isempty(new_particles)
        #= @show "start" =#
        particle_no += 1
        particle_inner = new_particles[1]
        particle_no % 200000 == 0 && @info "Particles Done: $(particle_no)"
        if length(new_particles) / og_num * model.n > 1000000000
            @info "Alot of ionisation events, consider stopping..."
            #= throw("Too many ionisation events") =#
        end
        #= event_number = 0 =#
        #= particle_states = fill(deepcopy(particle_inner), 3) =#
        #= dts_states = fill(fill(dt_EVENT(1.0 * s̃, :before, :before), 7), 3) =#
        # @show length(new_particles)
        while particle_inner.alive
            #= event_number += 1 =#
            #@debug "Particle before: " particle
            #' @debug debugParticle("before", particle_inner)
            #= checkBins(dt_EVENT(1.0 * s̃, :before, :before), particle_inner, params, particle_states, dts_states) =#
            checkBins(dt_EVENT(1.0 * s̃, :before, :before), particle_inner, params)
            dt = getDt!(dts, dt_times, particle_inner, model, params)

            :f in model.meas ? measureCumulative!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner) : nothing
            updateParticle!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner)
            #' @debug debugParticle("after flight", particle_inner)

            event!(dt, particle_inner, model, params, meas_inner, new_particles)
            #= particle_states[event_number%3+1] = deepcopy(particle_inner) =#
            #= dts_states[event_number%3+1] = dts =#

            #' @debug debugParticle("after event", particle_inner)
            #= @debug sum(meas_inner.measg.coll_count) > 5 && break =#

            particle_inner.alive = particle_inner.alive && checkAlive(model.type, particle_inner, model, meas_inner, params)
            measureEvent!(model.type, dt, particle_inner, model, params, meas_inner)
            # if verbose && meas_inner.measg.coll_count.elastic > 5 alive = false end
            #= checkBins(dt, particle_inner, params, particle_states, dts_states) =#
            checkBins(dt, particle_inner, params)
        end
        measureDeath!(model.type, particle_inner, model, params, meas_inner)
        popfirst!(new_particles)
    end
    return meas_inner::MEAS
end

function flightGetParticles(particle_og::Vector{PART}, model, params; coll_limit=Inf, t_limit=Inf * s)
    #' @debug "Conducting flight..."
    :f in model.meas && (params.a[1] > 0m^4 / s^2 || params.a[2] > 0m^4 / s^2) && @warn "Measuring flux with x or y field, third order will be wrong."
    meas_inner = initMeas(model, params)
    # Preallocate arrays
    #= dts = Array{dt_EVENT,1}(undef, 7) =#
    #= dt_times = Array{Quantity,1}(undef, 7) =#
    dts = repeat([dt_EVENT(1.0 * s̃, :before, :before)], 7)
    dt_times = repeat([1.0 * s̃], 7)

    measureTbin!(particle_og, model, params, meas_inner) ## measure initial stuff

    new_particles = particle_og
    i = 0
    while i != length(new_particles)
        i += 1
        particle_inner = new_particles[i]
        #= (length(new_particles)-i) % 1000 == 0 && @info length(new_particles) =#
        if length(new_particles) - i > 50000
            error("Too many ionisation events")
        end
        #= particle_states = fill(deepcopy(particle_inner), 3) =#
        #= dts_states = fill(fill(dt_EVENT(1.0 * s̃, :before, :before), 7), 3) =#
        # @show length(new_particles)
        while particle_inner.alive && sum(meas_inner.measg.coll_count) < coll_limit
            #@debug "Particle before: " particle
            #' @debug debugParticle("before", particle_inner)
            checkBins(dt_EVENT(1s, :before, :before), particle_inner, params)
            #= checkBins(dt_EVENT(1s, :before, :before), particle_inner, params, particle_states, dts_states) =#
            dt = getDt!(dts, dt_times, particle_inner, model, params)

            particle_inner.t + dt.time > model.type.t_f && @error debugParticle("Time error", particle_inner) dt
            :f in model.meas ? measureCumulative!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner) : nothing
            updateParticle!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner)

            event!(dt, particle_inner, model, params, meas_inner, new_particles)

            #' @debug debugParticle("after", particle_inner)
            #' @debug sum(meas_inner.measg.coll_count) > 5 && break

            particle_inner.alive = particle_inner.alive && checkAlive(model.type, particle_inner, model, meas_inner, params)
            measureEvent!(model.type, dt, particle_inner, model, params, meas_inner)
            # if verbose && meas_inner.measg.coll_count.elastic > 5 alive = false end
            checkBins(dt, particle_inner, params)
            #= checkBins(dt, particle_inner, params, particle_states, dts_states) =#
        end
        measureDeath!(model.type, particle_inner, model, params, meas_inner)
        #= popfirst!(new_particles) =#
    end
    return filter!((x) -> x.weight != 0.0, new_particles)
end


function flightGetParticles(particle_og::PART, model, params; coll_limit=Inf, t_limit=Inf * s)
    #' @debug "Conducting flight..."
    :f in model.meas && (params.a[1] > 0m^4 / s^2 || params.a[2] > 0m^4 / s^2) && @warn "Measuring flux with x or y field, third order will be wrong."
    meas_inner = initMeas(model, params)
    # Preallocate arrays
    #= dts = Array{dt_EVENT,1}(undef, 7) =#
    #= dt_times = Array{Quantity,1}(undef, 7) =#

    dts = repeat([dt_EVENT(1.0 * s̃, :before, :before)], 7)
    dt_times = repeat([1.0 * s̃], 7)
    measureTbin!(particle_og, model, params, meas_inner) ## measure initial stuff

    new_particles = [particle_og]
    i = 0
    while i != length(new_particles)
        i += 1
        particle_inner = new_particles[i]
        #= (length(new_particles)-i) % 1000 == 0 && @info length(new_particles) =#
        if length(new_particles) - i > 50000
            error("Too many ionisation events")
        end
        #= particle_states = fill(deepcopy(particle_inner), 3) =#
        #= dts_states = fill(fill(dt_EVENT(1.0 * s̃, :before, :before), 7), 3) =#
        # @show length(new_particles)
        while particle_inner.alive && sum(meas_inner.measg.coll_count) < coll_limit
            #@debug "Particle before: " particle
            #' @debug debugParticle("before", particle_inner)
            checkBins(dt_EVENT(1s, :before, :before), particle_inner, params)
            #= checkBins(dt_EVENT(1s, :before, :before), particle_inner, params, particle_states, dts_states) =#
            dt = getDt!(dts, dt_times, particle_inner, model, params)

            particle_inner.t + dt.time > model.type.t_f && @error debugParticle("Time error", particle_inner) dt
            :f in model.meas ? measureCumulative!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner) : nothing
            updateParticle!(dt.time, particle_inner, model.En0, model.n0, model, params, meas_inner)

            event!(dt, particle_inner, model, params, meas_inner, new_particles)

            #' @debug debugParticle("after", particle_inner)
            #' @debug sum(meas_inner.measg.coll_count) > 5 && break

            particle_inner.alive = particle_inner.alive && checkAlive(model.type, particle_inner, model, meas_inner, params)
            measureEvent!(model.type, dt, particle_inner, model, params, meas_inner)
            # if verbose && meas_inner.measg.coll_count.elastic > 5 alive = false end
            checkBins(dt, particle_inner, params)
            #= checkBins(dt, particle_inner, params, particle_states, dts_states) =#
        end
        measureDeath!(model.type, particle_inner, model, params, meas_inner)
        #= popfirst!(new_particles) =#
    end
    return filter!((x) -> x.weight != 0.0, new_particles)
end

function checkBins(dt, particle, params)
    # ebin is not a tracked quantity
    # if params.ebins != nothing && particle.ebin < length(params.ebins)
    #     (params.ebins[particle.ebin]-(eps()+1e-6)*eV) <= getEnergy(particle.v) <= (params.ebins[particle.ebin+1]+(eps()+1e-6)*eV) || error("""ebin error
    #         $(uconvert(eV,params.ebins[particle.ebin]))
    #         $(uconvert(eV,getEnergy(particle.v)))
    #         $(uconvert(eV,params.ebins[particle.ebin+1]))
    #         $(dt)
    #                                                                                                                                """)
    # end
    if params.cbins != nothing && particle.cbin < length(params.cbins)
        (params.cbins[particle.cbin] - (eps() + 1e-8) * eV) <= getEnergy(particle.v) <= (params.cbins[particle.cbin+1] + (eps() + 1e-8) * eV) || error("""CBin error
                        $(uconvert(eV,params.cbins[particle.cbin]))
                        $(uconvert(eV,getEnergy(particle.v)))
                        $(uconvert(eV,params.cbins[particle.cbin+1]))
                        $(dt)
                        """
        )
    end
    if params.zbins != nothing && 0 < particle.zbin < length(params.zbins)
        if (params.zbins[particle.zbin] - (1e-12 + eps()) * m̃) <= particle.r[3] <= (params.zbins[particle.zbin+1] + (1e-12 + eps()) * m̃) && particle.r[3] >= params.zbins[1]
        else
            error("""
                    zbin error
                  $(params.zbins[1])
                    $(params.zbins[particle.zbin] - (1e-12 + eps()) * m̃)
                    $(params.zbins[particle.zbin])
                    $(particle.r[3])
                    $(params.zbins[particle.zbin+1])
                    $(params.zbins[particle.zbin+1] + (1e-12 + eps()) * m̃)
                    $(dt)
                   """
            )
        end
    end
    if params.tbins != nothing && particle.tbin < length(params.tbins) && particle.t > params.tbins[1]
        (params.tbins[particle.tbin] - eps() * s) <= particle.t <= (params.tbins[particle.tbin+1] + eps() * s) || error("""tbin error
                    $(params.tbins[particle.tbin])
                    $(particle.t)
                    $(params.tbins[particle.tbin+1])
                    $(dt)
                    """)
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
    if params.tbins != nothing
        alive = alive && particle.t < model.type.t_f#*params.time_scale
        alive || @debug "Particle timed out"
        # alive = alive && sum(meas_inner.measg.coll_count) < 5
    end

    return alive && checkAliveModel(type, particle, model, meas_inner, params)
end

function checkAliveModel(type::SWARM, particle, model, meas_inner, params)
    return true
end

function checkAliveModel(type::LμJ, particle, model, meas_inner, params)
    if getEnergy(particle.v) < type.thresh#*params.energy_scale
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

    if getEnergy(particle.v) < type.thresh#*params.energy_scale
        @debug "Absorbed"
    end

    if params.zbins != nothing
        # :SP_z in model.meas ? alive = alive && particle.zbin != length(params.zbins) : alive = alive && particle.zbin != length(params.zbins) && particle.zbin != 0
        #= alive = alive && particle.zbin != length(params.zbins) && particle.zbin > 1 =#
        alive = alive && particle.zbin != length(params.zbins) && particle.zbin != 1
        (particle.zbin == length(params.zbins) || particle.zbin == 1) && @debug "Detected Z"
    end

    if params.ybins != nothing
        alive = alive && particle.ybin != length(params.ybins)
        particle.ybin == length(params.ybins) && @debug "Detected Y"
    end

    if params.xbins != nothing
        alive = alive && particle.xbin != length(params.xbins)
        particle.xbin == length(params.xbins) && @debug "Detected X"
    end

    :SP in model.meas ? alive = alive && sum(meas_inner.measg.coll_count) == 0 : nothing

    return alive && getEnergy(particle.v) > type.thresh#*params.energy_scale
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
    if dt.type == :dt_zbin
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
    meas_inner.measc.r[particle.tbin, :] += particle.weight * (particle.r + 0.5 * particle.v * dt + (1 / 6) * a * (dt^2)) * dt
    meas_inner.measc.r²[particle.tbin, :] += @. particle.weight * dt * (
                                                    a .^ 2 * dt .^ 4 / 20 +
                                                    dt^2 * (a .* particle.r + particle.v .^ 2) / 3 +
                                                    a * dt .^ 3 .* particle.v / 4 +
                                                    particle.r .^ 2 +
                                                    particle.r .* particle.v * dt
                                                )

    meas_inner.measc.v[particle.tbin, :] += particle.weight * (particle.v + 0.5 * a * dt) * dt
    meas_inner.measc.rv[particle.tbin, :] += particle.weight * (particle.r .* particle.v + 0.5 * (particle.r .* a + particle.v .^ 2) * dt + 0.5 * particle.v .* a * (dt^2) + (1 / 8) * ((a) .^ 2) * (dt^3)) * dt
    # NOTE: These have been optimised for electric field in z direction
    #
    meas_inner.measc.r²v[particle.tbin, :] += 1 / 24 * dt * particle.weight * (
                                                  a * dt + 2 * particle.v) .* (
        6 * particle.r * dt .* (a * dt + 2 * particle.v) + dt^2 * (a * dt + 2particle.v) .^ 2 + 12particle.r .^ 2)

    meas_inner.measc.rx²vz[particle.tbin] += 1 / 12 * dt^2 * particle.weight * (
                                                 4 * dt * particle.v[1] * (2 * a[3] * particle.r[1] + particle.v[3] * particle.v[1]) +
                                                 6 * particle.r[1] * (a[3] * particle.r[1] + 2 * particle.v[1] * particle.v[3]) +
                                                 3 * a[3] * dt^2 * particle.v[1]^2 +
                                                 12 * particle.v[3] * particle.r[1]^2 / dt)

    meas_inner.measc.rzrxvx[particle.tbin] += particle.weight * particle.v[1] * dt * (
                                                  a[3] * particle.v[1] * dt^3 / 8 +
                                                  dt^2 * (a[3] * particle.r[1] + 2 * particle.v[1] * particle.v[3]) / 6 +
                                                  dt * (particle.v[1] * particle.r[3] + particle.v[3] * particle.r[1]) / 2 +
                                                  particle.r[3] * particle.r[1])

    meas_inner.measc.weight[particle.tbin] += particle.weight * dt
    meas_inner.measc.dt[particle.tbin] += dt
end

function measureDeath!(type::SWARM, particle, model, params, meas_inner)
end

function measureDeath!(type::LμJ, particle, model, params, meas_inner)
    if particle.rbin == length(params.rbins) && sum(meas_inner.measg.coll_count) != 0
        :χ in model.meas ? meas_inner.measLμJ.χ[findMinInd(params.χbins, abs(atan(rad, particle.r[2], particle.r[3])))] += particle.weight : nothing

        if :εχ in model.meas
            @debug """
            Velocity angle:  $(abs(atan(°,particle.v[2],particle.v[3])))
            Spatial angle:  $(abs(atan(°,particle.r[2],particle.r[3])))
            Velocity angle:  $(atan(°,particle.v[2],particle.v[3]))
            Spatial angle:  $(atan(°,particle.r[2],particle.r[3]))
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
            Velocity angle:  $(abs(atan(°,particle.v[2],particle.v[3])))
            Spatial angle:  $(abs(atan(°,particle.r[2],particle.r[3])))
            Velocity angle:  $(atan(°,particle.v[2],particle.v[3]))
            Spatial angle:  $(atan(°,particle.r[2],particle.r[3]))
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
    if :tz in model.meas && particle.zbin == length(params.zbins)
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
    if :SP in model.meas && sum(meas_inner.measg.coll_count) == 1
        meas_inner.measg.SP_ε = energy_loss
        meas_inner.measg.SP_r = sqrt(dot(particle.r, particle.r)) # Assumes particle started at [0,0,0]
    end
end

function event!(dt, particle, model, params, meas_inner, new_particles)
    if dt.type == :dt_collision
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
    measureTbin!(particle, model, params, meas_inner)
end

function measureTbin!(particle, model, params, meas_inner)
    if :t ∈ model.meas
        @. meas_inner.meast.v²[particle.tbin, :] += particle.v^2 * particle.weight
        @. meas_inner.meast.r²[particle.tbin, :] += particle.r^2 * particle.weight
        @. meas_inner.meast.r³[particle.tbin, :] += particle.r^3 * particle.weight
        @. meas_inner.meast.r[particle.tbin, :] += particle.r * particle.weight
        @. meas_inner.meast.v[particle.tbin, :] += particle.v * particle.weight
        @. meas_inner.meast.rv[particle.tbin, :] += particle.v * particle.r * particle.weight
        meas_inner.meast.zx²[particle.tbin] += particle.r[3] * particle.r[2]^2 * particle.weight
        meas_inner.meast.w[particle.tbin] += particle.weight
        meas_inner.meast.εᵢ[particle.tbin] += 0.5 * me * dot(particle.v, particle.v) * particle.weight
    end

    :zyt ∈ model.meas ? meas_inner.measLμJ.zyt[particle.zbin, particle.ybin, particle.tbin] += particle.weight : nothing

    :tz ∈ model.meas ? meas_inner.measCell.tz[particle.tbin, particle.zbin] += particle.weight : nothing
    :tz_v ∈ model.meas ? meas_inner.measCell.v_tz[particle.tbin, particle.zbin] += particle.weight * particle.v[3] : nothing
    # :tz ∈ model.meas && params.zbins[2] < particle.r[3] < params.zbins[end] ? meas_inner.measCell.tz[particle.tbin,1] += particle.weight : nothing
    # :tz ∈ model.meas && particle.r[3] < params.zbins[end] ? meas_inner.measCell.tz[particle.tbin,1] += particle.weight : nothing
end


#### * Collisions

function eventCollision!(dt::dt_EVENT, particle, model, params, meas_inner, new_particles)
    new_particle = nothing
    cs_ind = getCollType(particle, model, params)
    energy_loss = 0.0
    :SP in model.meas && sum(meas_inner.measg.coll_count) == 0 ? energy_loss = getEnergy(particle.v) : nothing
    #' @debug """
    #' Collision: 
    #' $(cs_ind == 0 ? "NULL" : params.gas[cs_ind][1].type) $(dt.time)
    #' """
    if cs_ind != 0
        cs, weight, mass = params.gas[cs_ind]
        g = getRelVel(cs_ind, cs, mass, particle, model, params)

        if cs.type == :ELASTIC
            elasticColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
        elseif cs.type == :EXCITATION
            excitationColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
            # energy_loss = cs.energy_loss
        elseif cs.type == :IONIZATION
            ionisationColl!(cs_ind, cs, mass, model.ion_type, g, particle, model, params, meas_inner, new_particles)
            # energy_loss = cs.energy_loss
        elseif cs.type == :ATTACHMENT
            attachmentColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
        end
    end

    :SP in model.meas && sum(meas_inner.measg.coll_count) == 1 ? energy_loss -= getEnergy(particle.v) : nothing

    measureCollision!(energy_loss, particle, model, params, meas_inner)
    updateCbin!(particle, params)
    return new_particle
end

function getRelVel(cs_no, cs, mass, particle, model, params)
    if model.T.val == 0.0 * K
        g = copy(particle.v)
    else
        g = findu(cs_no, cs, mass, particle, model, params)
    end
    #' @debug g
    return g
end

function getCollType(particle, model, params)
    begin
        R = rand()
        CS = 0.0Å^2 * sqrt(eV)
        en = uconvert(eV, getEnergy(particle.v))#/params.energy_scale
        #' v = sqrt(dot(particle.v, particle.v))#/params.length_scale

        for (cs_ind, (cs, weight, mass)) in enumerate(params.gas)
            CS += CollisionFrequency(cs)(sqrt(en)) * weight

            if R < CS / getTotalFreq(particle, model, params)
                return cs_ind
            end
        end
        return 0
    end
end

function elasticColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)
    #me2 = 4amu
    # @debug uconvert(amu,me), uconvert(amu,mass)
    coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
    v2 = norm(g) * SArray{Tuple{3},Float64,1,3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
    v2 = shift(v2, g)
    v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)
    meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
    meas_inner.measg.coll_count[cs_no] += 1
    # meas_inner.no_elas += 1
    particle.v = v2
end

function excitationColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)

    if dot(g, g) - 2 * cs.energy_loss / me > 0m^2 / s^2
        coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
        v2 = sqrt(dot(g, g) - 2 * cs.energy_loss / me) * SArray{Tuple{3},Float64,1,3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
        v2 = shift(v2, g)
        v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)
        meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
        meas_inner.measg.coll_count[cs_no] += 1
        particle.v = v2
    end

end

function ionisationColl!(cs_no, cs, mass, ion_type, g, particle, model, params, meas_inner, new_particles)

    if dot(g, g) - 2 * cs.energy_loss / me > 0m^2 / s^2
        part_energy, other_part_energy = getIonisationEnergy(cs.energy_loss, ion_type, g, particle, model, params)

        coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)

        v2 = sqrt(2 * part_energy / me) * SArray{Tuple{3},Float64,1,3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
        v2 = shift(v2, g)
        v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)

        meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
        meas_inner.measg.coll_count[cs_no] += 1
        particle.v = v2


        if ion_type.method == :WEIGHTING
            particle.type == :electron ? particle.weight *= 2.0 : nothing
        elseif ion_type.method == :EXPLICIT
            coschi2, sinchi2, phi2 = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
            v2 = sqrt(2 * other_part_energy / me) * SArray{Tuple{3},Float64,1,3}(cos(phi2) * sinchi2, sin(phi2) * sinchi2, coschi2)
            new_particle = deepcopy(particle)
            meas_inner.measg.coschi += dot(new_particle.v, v2) / (norm(v2) * norm(new_particle.v))
            meas_inner.measg.coll_count[cs_no] += 1
            new_particle.v = v2
            # @info debugParticle("other particle",new_particle)
            updateCbin!(new_particle, params)
            push!(new_particles, new_particle)
        else
            error("Unknown non-conservative method in ion_type")
        end
    else
        @error begin
            """
            Ion threshold greater than energy"
            Cross section #: $(cs_no)
            Threshold: $(cs.energy_loss)
            Dot: $(dot(g,g))
            minus: $(2*cs.energy_loss/me)
            $(debugParticle("Ionisation collision",particle))

            """
        end
    end

end

function getIonisationEnergy(thresh, ion_type::DiscreteIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type == :positron
        other_part_energy = ion_type.loss#*params.energy_scale*q
    else
        other_part_energy = total_energy - ion_type.loss#*params.energy_scale*q
    end
    part_energy = total_energy - other_part_energy
    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::UniformIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    other_part_energy = total_energy * rand()
    part_energy = total_energy - other_part_energy
    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::RatioIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type == :positron
        other_part_energy = total_energy * ion_type.ratio(getEnergy(g))::Float64
        # other_part_energy = 0
    else
        other_part_energy = total_energy * (1 - ion_type.ratio(getEnergy(g))::Float64)
    end
    part_energy = total_energy - other_part_energy
    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::FunctionalIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type == :positron
        other_part_energy = ion_type.func[particle.cbin](rand())#*params.energy_scale*q
    else
        other_part_energy = total_energy - ion_type.func[particle.cbin](rand())#*params.energy_scale*q
    end
    part_energy = total_energy - other_part_energy

    if part_energy < 0
        @info "part" part_energy#/params.energy_scale/q
        @info "other" other_part_energy#/params.energy_scale/q
    end

    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::DiscreteFunctionalIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type == :positron
        other_part_energy = ion_type.loss(getEnergy(g))#/params.energy_scale)*params.energy_scale*q
    else
        other_part_energy = total_energy - ion_type.loss(getEnergy(g))#/params.energy_scale)*params.energy_scale*q
    end
    part_energy = total_energy - other_part_energy

    if part_energy < 0
        @info "part" part_energy#/params.energy_scale/q
        @info "other" other_part_energy#/params.energy_scale/q
    end

    return part_energy, other_part_energy
end

function attachmentColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)
    if model.ion_type.method == :WEIGHTING
        particle.type == :positron ? particle.alive = false : particle.weight /= 2.0
    elseif model.ion_type.method == :EXPLICIT
        particle.alive = false
        particle.weight = 0.0
    else
        error("Unknown non-conservative method in ion_type")
    end
end

function scatteringAngle(cs_no, cs, mass, g, particle, model, params)
    coschi = 1 - 2 * rand()
    #' coschi = params.splat_gas.dcs[cs_no][particle.cbin](rand())::Float64
    #coschi = -1
    sinchi = sqrt(1 - coschi^2)
    phi = rand() * 2 * pi
    #' @debug """
    #' Scattering angle: 
    #' $coschi $sinchi $phi"
    #' """
    return coschi, sinchi, phi
end

#@fastmath function shift(v_new::SVector{3, Quantity},v::SVector{3, Quantity})
@fastmath function shift(v_new, v)
    relvel_mag = norm(v)
    if relvel_mag == 0
        orig_costheta = 1.0
        orig_sintheta = 0.0
        orig_cosphi = 1.0
        orig_sinphi = 0.0
    else
        orig_costheta = v[3] / relvel_mag
        orig_sintheta = sqrt(1 - orig_costheta^2)
        if orig_sintheta == 0
            orig_cosphi = 1.0
            orig_sinphi = 0.0
        else
            orig_cosphi = v[1] / relvel_mag / orig_sintheta
            orig_sinphi = v[2] / relvel_mag / orig_sintheta
        end
    end

    Ry = @SMatrix [orig_costheta 0.0 orig_sintheta; 0.0 1.0 0.0; -orig_sintheta 0.0 orig_costheta]
    Rz = @SMatrix [orig_cosphi -orig_sinphi 0.0; orig_sinphi orig_cosphi 0.0; 0.0 0.0 1.0]

    return Rz * Ry * v_new
end
#### * Update Bin numbers

function updateCbin!(particle, params)
    ea = getEnergy(particle.v)
    while ea < params.cbins[particle.cbin]
        particle.cbin -= 1
    end
    if particle.cbin != length(params.cbins)
        while ea > params.cbins[particle.cbin+1]
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

function eventEnergy!(dt::dt_EVENT, particle, model, params, meas_inner)
    if dt.direction == :Up
        #' @debug "Energy Bin Up: #$(particle.cbin) $(dt.time)"
        particle.cbin += 1
    elseif dt.direction == :Down
        #' @debug "Energy Bin Down: #$(particle.cbin) $(dt.time)"
        particle.cbin -= 1
    end
end

function eventRBin!(dt::dt_EVENT, particle, model, params, meas_inner)
    if dt.direction == :Up
        #' @debug "RBin Up: #$(particle.rbin) $(dt.time)"
        if particle.r[2]^2 + particle.r[3]^2 < params.rbins[particle.rbin+1]
            # If particle is going up an Rbin, it can't be less than that bin value
        end
        particle.rbin += 1
    elseif dt.direction == :Down
        #' @debug "RBin Down: #$(particle.rbin) $(dt.time)"
        if particle.r[2]^2 + particle.r[3]^2 > params.rbins[particle.rbin]
            # If particle is going down an Rbin, it can't be more than the existing bin value
        end
        particle.rbin -= 1
    end
end

function eventYBin!(dt::dt_EVENT, particle, model, params, meas_inner)
    if dt.direction == :Up
        #' @debug "YBin Down: #$(particle.ybin) $(dt.time)"
        particle.ybin += 1
    elseif dt.direction == :Down
        #' @debug "YBin Down: #$(particle.ybin) $(dt.time)"
        particle.ybin -= 1
    end
end

function eventZBin!(dt::dt_EVENT, particle, model, params, meas_inner)
    if dt.direction == :Up
        #' @debug "ZBin Up: #$(particle.zbin) $(dt.time)"
        particle.zbin += 1
    elseif dt.direction == :Down
        #' @debug "ZBin Down: #$(particle.zbin) $(dt.time)"
        particle.zbin -= 1
    end
    if :SP_z in model.meas && particle.zbin != 0
        meas_inner.measg.SP_εs[particle.zbin] = getEnergy(particle.v)
    end
end

function eventXBin!(dt::dt_EVENT, particle, model, params, meas_inner)
    if dt.direction == :Up
        #' @debug "XBin Down: #$(particle.xbin) $(dt.time)"
        particle.xbin += 1
    elseif dt.direction == :Down
        #' @debug "XBin Down: #$(particle.xbin) $(dt.time)"
        particle.xbin -= 1
    end
end

# * Fetch Functions

#### * Temperature
@fastmath function findu(cs_ind, cs, mass, particle, model, params)
    #TODO: make generic function based on inputs
    w = sqrt(2 * k * model.T.val / mass)
    v = particle.v#*params.time_scale/params.length_scale
    he = electron_v2eV * (norm(v) + 5.0 * w)^2
    le = max(1e-25u"kg*m^2/s^2", electron_v2eV * (norm(v) - 5.0 * w) .^ 2)

    #TODO: get actual CF from CS

    hcf = evalTotalFreq(params.gas, he)
    lcf = evalTotalFreq(params.gas, params)

    g = copy(v)
    mcf = max(hcf, lcf)
    while true
        # usually it only goes through one while loop iteration
        u = SArray{Tuple{3},Float64,1,3}(randn(3)) * w ./ sqrt(2.0) # could only run this
        g = v .- u
        # @show v
        # @show u
        # @show g
        eg = electron_v2eV * dot(g, g)
        # @show eg
        ccf = evalTotalFreq(params.gas, eg)
        # @show lcf
        # @show hcf
        # @show ccf
        # @show mcf

        if ccf / mcf == 0.0
            throw("Zero relative energy")
        end
        if rand() < ccf / mcf
            break
        end
    end
    #' @debug "Neutral velocity: x = $(g[1]), y = $(g[2]), z = $(g[3])"
    return g - +SArray{Tuple{3},Float64,1,3}(1e2, 0, 0)m / s#*params.length_scale/params.time_scale
end
