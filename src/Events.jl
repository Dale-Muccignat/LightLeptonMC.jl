function eventEnergy!(dt::dt_EVENT, particle, model, params, meas_inner)
    return if dt.direction == :Up
        #' @debug "Energy Bin Up: #$(particle.cbin) $(dt.time)"
        particle.cbin += 1
    elseif dt.direction == :Down
        #' @debug "Energy Bin Down: #$(particle.cbin) $(dt.time)"
        particle.cbin -= 1
    end
end

function eventRBin!(dt::dt_EVENT, particle, model, params, meas_inner)
    return if dt.direction == :Up
        #' @debug "RBin Up: #$(particle.rbin) $(dt.time)"
        if particle.r[2]^2 + particle.r[3]^2 < params.rbins[particle.rbin + 1]
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
    return if dt.direction == :Up
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
    return if :SP_z in model.meas && particle.zbin != 0
        meas_inner.measg.SP_εs[particle.zbin] = getEnergy(particle.v)
    end
end

function eventXBin!(dt::dt_EVENT, particle, model, params, meas_inner)
    return if dt.direction == :Up
        #' @debug "XBin Down: #$(particle.xbin) $(dt.time)"
        particle.xbin += 1
    elseif dt.direction == :Down
        #' @debug "XBin Down: #$(particle.xbin) $(dt.time)"
        particle.xbin -= 1
    end
end
