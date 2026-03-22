#TODO: Add TimerOutputs
#TODO: Only select electrons with position and velocity directed toward detector in LuJ

abstract type Direction end
abstract type Event end
struct Down <: Direction end
struct Up <: Direction end

struct dt_EVENT
    time::uT
    type::Symbol
    direction::Symbol
end

function getDt(particle, model, params)
    min_dt = dt_EVENT(typemax(uT), :NONE, :NONE)

    @inbounds begin
        candidates = (
            getCDt(particle, model, params),
            getEnDt(model.En0, model.n0, particle, model, params),
            getTDt(particle, model, params),
            getRDt(particle, model, model.type, params),
            getZDt(model.En0, particle, model, params),
            getYDt(particle, model, params),
            getXDt(particle, model, params),
        )

        for dt in eachindex(candidates)
            if candidates[dt].time < min_dt.time
                min_dt = candidates[dt]
            end
        end
    end

    return min_dt
end

function getEnDt(En0::UniformField, n0::UniformDensity, particle, model, params)
    if params.cbins |> !isempty
        v_perp = particle.v[1]^2 + particle.v[2]^2 #Seperate out of the vz velocity
        if (particle.v[3] < 0m / s) && (((params.cbins[particle.cbin]) / electron_v2eV - v_perp) > 0 * (m / s)^2) && (particle.cbin != 1)
            # if (particle.v[3] < 0m/s) && (((params.cbins[particle.cbin]-1e-6eV)/electron_v2eV - v_perp) > 0*(m/s)^2) && (particle.cbin!= 1)
            #If the velocity is less than zero, i.e if the energy bin needs to go down
            #The second condition is to ensure that the energy can go down considering the perpindicular velocity... or something along those lines. It just makes sure the squareroot doesn't return an imaginary number
            vf_z = sqrt((params.cbins[particle.cbin]) / electron_v2eV - v_perp)
            # vf_z is the velocity the particle will need to reach the lower boundry of the energy bin
            return dt_EVENT(-(particle.v[3] + vf_z) / (params.a[3] * params.n0), :dt_energy, :Down)
            #Just finding the time until the particles velocity reaches vf_z
        else #if (particle.cbin_no < length(vars.cbins))
            # here the vz is positive, and also making sure we aren't at the max bin value
            #Don't need to check the second  condition above as the particle can only increase in energy
            vf_z = sqrt((params.cbins[particle.cbin + 1]) / electron_v2eV - v_perp)
            #same as before except with the lower boundry of the next bin up
            return dt_EVENT(-(particle.v[3] - vf_z) / (params.a[3] * params.n0), :dt_energy, :Up)
        end
    end
    return dt_EVENT(Inf * s, :dt_energy, :Up)
end

function getEnDt(Eno::ZeroField, n0, particle, model, params)
    return dt_EVENT(Inf * s, :dt_energy, :Up)
end

@fastmath function getRDt(particle, model, type::LμJ, params)
    #NOTE: Only for fields of 0
    up, down = Inf * s, Inf * s
    if params.rbins !== nothing
        # Equation of circle
        # r²=x²+y²
        # Projectile motion
        # rₓ=rₓ₀+vₓ₀t
        # rₓ=rₓ₀+vₓ₀t
        #
        p = particle.r
        # if particle.rbin !== 1
        #if particle is inside detector of jet, the time to bin down should be offset by jet_center
        p -= model.type.jet_center
        # end
        v = particle.v
        r21 = params.rbins[particle.rbin]
        r22 = params.rbins[particle.rbin + 1]
        a = v[3]^2 + v[2]^2
        b = 2 * (p[3] * v[3] + v[2] * p[2])
        c = p[3]^2 + p[2]^2 - r21 + 1.0e-16m^2
        delta = b^2 - 4 * a * c
        t1, t2 = Inf * s, Inf * s
        if (delta > 0m^4 / s^2) #&& (c > 0.)
            t1 = (-b + sqrt(delta)) / (2 * a)
            t2 = (-b - sqrt(delta)) / (2 * a)
            if t1 <= 0s
                t1 = Inf * s
            end
            if t2 <= 0s
                t2 = Inf * s
            end
            down = min(t1, t2)
            if r21 == 0s
                down = Inf * s
            end
        end
        @debug """
        RBin Debug
        shifted particle: $p
        velocity: $v
        rs: $r21 $r22
        ts: $t1 $t2
        """
        # if particle.rbin !== 1
        # undo offset from above if next moving up to detector
        # p += model.type.jet_center
        # end
        if particle.rbin == length(params.rbins) - 1
            # offset if moving up to jet edge
            p += model.type.jet_center
        end
        b = 2 * (p[3] * v[3] + v[2] * p[2])
        c = p[3]^2 + p[2]^2 - r22 - 1.0e-16m^2
        delta = b^2 - 4 * a * c
        if (delta > 0m^4 / s^2) #&& (c < 0.)
            t1 = (-b + sqrt(delta)) / (2 * a)
            t2 = (-b - sqrt(delta)) / (2 * a)
            if t1 <= 0s
                t1 = Inf * s
            end
            if t2 <= 0s
                t2 = Inf * s
            end
            up = min(t1, t2)
        end
        @debug """
        RBin Debug
        shifted particle: $p
        velocity: $v
        rs: $r21 $r22
        ts: $t1 $t2
        """
    end
    if up == min(up, down)
        return dt_EVENT(up, :dt_rbin, :Up)
    else
        return dt_EVENT(down, :dt_rbin, :Down)
    end
end
@fastmath function getRDt(particle, model, type, params)
    #NOTE: Only for fields of 0
    up, down = Inf * s, Inf * s
    if params.rbins |> !isempty
        p = particle.r
        v = particle.v
        r21 = params.rbins[particle.rbin]
        r22 = params.rbins[particle.rbin + 1]
        a = v[3]^2 + v[2]^2
        b = 2 * (p[3] * v[3] + v[2] * p[2])
        c = p[3]^2 + p[2]^2 - r21 + 1.0e-16m^2
        delta = b^2 - 4 * a * c
        if (delta > 0m^4 / s^2) #&& (c > 0.)
            t1 = (-b + sqrt(delta)) / (2 * a)
            t2 = (-b - sqrt(delta)) / (2 * a)
            if t1 <= 0s
                t1 = Inf * s
            end
            if t2 <= 0s
                t2 = Inf * s
            end
            down = min(t1, t2)
            if r21 == 0s
                down = Inf * s
            end
        end
        c = c + r21 - r22 - 2.0e-16m^2
        delta = b^2 - 4 * a * c
        if (delta > 0m^4 / s^2) #&& (c < 0.)
            t1 = (-b + sqrt(delta)) / (2 * a)
            t2 = (-b - sqrt(delta)) / (2 * a)
            if t1 <= 0s
                t1 = Inf * s
            end
            if t2 <= 0s
                t2 = Inf * s
            end
            up = min(t1, t2)
        end
    end
    if up == min(up, down)
        return dt_EVENT(up, :dt_rbin, :Up)
    else
        return dt_EVENT(down, :dt_rbin, :Down)
    end
end

function getCDt(particle, model, params)
    return dt_EVENT(randexp() / (getTotalFreq(particle, model, params) / sqrt(electron_v2eV) * getDensity(particle, model.n0, params)), :dt_collision, :Up)
end

function getZDt(field::ZeroField, particle, model, params)
    if params.zbins |> !isempty
        if model.optimise_spatial == true
            if (particle.v[3] < 0m̃ / s̃) && (particle.zbin != 0)
                return dt_EVENT((params.zbins[1] - particle.r[3]) / particle.v[3], :dt_zbin, :Down)
            elseif (particle.zbin < length(params.zbins)) && (particle.v[3] > 0m / s)
                return dt_EVENT((params.zbins[end] + 1.0e-12m̃ - particle.r[3]) / particle.v[3], :dt_zbin, :Up)
            end
        else
            if (particle.v[3] < 0m̃ / s̃) && (particle.zbin != 0)
                return dt_EVENT((params.zbins[particle.zbin] - 1.0e-12m̃ - particle.r[3]) / particle.v[3], :dt_zbin, :Down)
            elseif (particle.zbin < length(params.zbins)) && (particle.v[3] > 0m / s)
                return dt_EVENT((params.zbins[particle.zbin + 1] + 1.0e-12m̃ - particle.r[3]) / particle.v[3], :dt_zbin, :Up)
            end
        end
    end
    return dt_EVENT(Inf * s, :dt_zbin, :Up)
end

function getZDt(field::UniformField, particle, model, params)
    if params.zbins |> !isempty
        # if (particle.v[3] < 0m/s) && (( (particle.v[3]^2)/(params.a[3]*params.n0) + 2*params.zbins[particle.zbin]-2e-10m ) > 2*particle.r[3]) && (particle.zbin != 0)
        if (particle.v[3] < 0m / s) && ((particle.v[3])^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin] + 1.0e-12m̃) > 0m^2 / s^2) && (particle.zbin != 0)
            # return dt_EVENT((- particle.v[3] - sqrt( (particle.v[3])^2 - 2*(params.a[3]*params.n0)*(particle.r[3]-params.zbins[particle.zbin]+1e-11m) ))/(params.a[3]*params.n0), :dt_zbin,:Down)
            a = (-particle.v[3] - sqrt((particle.v[3])^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin] + 1.0e-12m̃))) / (params.a[3] * params.n0)
            b = (-particle.v[3] + sqrt((particle.v[3])^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin] + 1.0e-12m̃))) / (params.a[3] * params.n0)
            b > 0s && b < a && error("Zbin calculation error down: a=$(a), b=$(b)")
            # if particle.v[3]^2 - 2*(params.a[3]*params.n0)*(particle.r[3]-params.zbins[particle.zbin+1]-1e-14m̃ ) > 0m̃^2/s̃^2
            #     c = (- particle.v[3] + sqrt( particle.v[3]^2 - 2*(params.a[3]*params.n0)*(particle.r[3]-params.zbins[particle.zbin+1]-1e-14m̃ ) ))/(params.a[3]*params.n0)
            #     d = (- particle.v[3] - sqrt( particle.v[3]^2 - 2*(params.a[3]*params.n0)*(particle.r[3]-params.zbins[particle.zbin+1]-1e-14m̃ ) ))/(params.a[3]*params.n0)
            #     (c > 0s && c < a) || (d > 0s && d < a) && error("Zbin calculation error down: a=$(a), b=$(b), c=$(c), d=$(d)")
            # end
            return dt_EVENT(a, :dt_zbin, :Down)
        elseif (particle.zbin < length(params.zbins))
            # @debug particle.r[3] particle.v[3] params.zbins[particle.zbin:(particle.zbin+1)]
            # return dt_EVENT((- particle.v[3] + sqrt( particle.v[3]^2 - 2*(params.a[3]*params.n0)*(particle.r[3]-params.zbins[particle.zbin+1]-1e-10m) ))/(params.a[3]*params.n0), :dt_zbin, :Up)
            if (particle.v[3]^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin + 1])) < 0m^2 / s^2
                @info debugParticle("Zbin error", particle)
                @info params.zbins[particle.zbin]
                @info params.zbins[particle.zbin + 1]
            end
            a = (-particle.v[3] + sqrt(particle.v[3]^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin + 1] - 1.0e-12m̃))) / (params.a[3] * params.n0)
            b = (-particle.v[3] - sqrt(particle.v[3]^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin + 1] - 1.0e-12m̃))) / (params.a[3] * params.n0)
            b > 0s && b < a && error("Zbin calculation error up: a=$(a), b=$(b)")
            if (particle.v[3])^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin] + 1.0e-11m̃) > 0m̃^2 / s̃^2
                c = (-particle.v[3] - sqrt((particle.v[3])^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin] + 1.0e-12m̃))) / (params.a[3] * params.n0)
                d = (-particle.v[3] + sqrt((particle.v[3])^2 - 2 * (params.a[3] * params.n0) * (particle.r[3] - params.zbins[particle.zbin] + 1.0e-12m̃))) / (params.a[3] * params.n0)
                (c > 0s && c < a) || (d > 0s && d < a) && error(
                    """ Zbin calculation error up: a=$(a), b=$(b), c=$(c), d=$(d)
                    $(debugParticle("Zbin error", particle))
                    $(params.zbins[particle.zbin])
                    $(params.zbins[particle.zbin + 1])
                    """
                )
            end
            return dt_EVENT(a, :dt_zbin, :Up)
        end
    end
    return dt_EVENT(Inf * s, :dt_zbin, :Up)
end

function getYDt(particle, model, params)
    if params.ybins |> !isempty
        if model.optimise_spatial == true
            if (particle.v[2] < 0m / s) && (particle.ybin != 0)
                return dt_EVENT((params.ybins[1] - particle.r[2]) / particle.v[2], :dt_ybin, :Down)
            elseif (particle.ybin < length(params.ybins)) && (particle.v[2] > 0m / s)
                return dt_EVENT((params.ybins[end] + 1.0e-12m - particle.r[2]) / particle.v[2], :dt_ybin, :Up)
            end
        else
            if params.ybins |> !isempty
                if (particle.v[2] < 0m / s) && (particle.ybin != 0)
                    return dt_EVENT((params.ybins[particle.ybin] - 1.0e-12m - particle.r[2]) / particle.v[2], :dt_ybin, :Down)
                elseif (particle.ybin < length(params.ybins)) && (particle.v[2] > 0.0)
                    return dt_EVENT((params.ybins[particle.ybin + 1] + 1.0e-12m - particle.r[2]) / particle.v[2], :dt_ybin, :Up)
                end
            end
        end
    end
    return dt_EVENT(Inf * s, :dt_ybin, :Up)
end

function getXDt(particle, model, params)
    if params.xbins |> !isempty
        if model.optimise_spatial == true
            if (particle.v[1] < 0m / s) && (particle.xbin != 0)
                return dt_EVENT((params.xbins[1] - particle.r[1]) / particle.v[1], :dt_xbin, :Down)
            elseif (particle.xbin < length(params.xbins)) && (particle.v[1] > 0.0)
                return dt_EVENT((params.xbins[end] + 1.0e-12m - particle.r[1]) / particle.v[1], :dt_xbin, :Up)
            end
        else
            if params.xbins |> !isempty
                if (particle.v[1] < 0m / s) && (particle.xbin != 0)
                    return dt_EVENT((params.xbins[particle.xbin]m - particle.r[1]) / particle.v[1], :dt_xbin, :Down)
                elseif (particle.xbin < length(params.xbins)) && (particle.v[1] > 0.0)
                    return dt_EVENT((params.xbins[particle.xbin + 1] + 1.0e-12m - particle.r[1]) / particle.v[1], :dt_xbin, :Up)
                end
            end
        end
    end
    return dt_EVENT(Inf * s, :dt_xbin, :Up)
end

function getTDt(particle, model, params)
    if params.tbins |> !isempty
        if particle.tbin < length(params.tbins) && particle.t >= params.tbins[1]
            return dt_EVENT(params.tbins[particle.tbin + 1] - particle.t, :dt_tbin, :Up)
        end
    end

    return dt_EVENT(Inf * s̃, :dt_tbin, :Up)

end
