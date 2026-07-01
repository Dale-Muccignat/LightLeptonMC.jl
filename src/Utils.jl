module Utils

using ..Structs
using MCSimUnits, LinearAlgebra, Roots, AndExport
export debugDts, debugParticle, getEnergy, getVelocity, structprint, adaptiveMesh, adaptiveMeshlog

const electron_v2eV = 0.5 * me

@xport function getTotalFreq(particle, model, params)
    params.total_freq[particle.cbin]
end

@xport function getDensity(particle, n0::UniformDensity, params)
    params.n0
end
@xport function getDensity(particle, n0::LμJDensity, params)
    params.n0[particle.rbin]
end
@xport function getDensity(particle, n0::LμJDensityGradient, params)
    params.n0[particle.rbin]
end
@xport function getDensity(particle, n0::CellDensity, params)
    params.n0[particle.zbin]
end

@xport function findTotalCF(en, model, params)
    CS = 0.0 * u"Å^2*m/s"
    for cs_ind in eachindex(params.splat_gas.cs)
        @inbounds CS += params.splat_gas.cs[cs_ind](en) * params.splat_gas.w[cs_ind] * sqrt(en / electron_v2eV)
    end

    return CS
end

function debugParticle(header, particle)
    """
    Particle $header: 
      ε: $(uconvert(eV,getEnergy(particle.v)))
      t: $(particle.t)
    pos: x=$(particle.r[1]) y=$(particle.r[2]) z=$(particle.r[3]) 
    vel: x=$(particle.v[1]) y=$(particle.v[2]) z=$(particle.v[3]) 
    rad: $(sqrt(particle.r[2]^2+particle.r[3]^2))
    weight: $(particle.weight) alive: $(particle.alive)
    bin: t=$(particle.tbin) | c=$(particle.cbin) | e=$(particle.ebin) | z=$(particle.zbin) 
    """
end

function debugDts(dts)
    """ time of flights
        CS Bin: $(dts[1].time)
        En Bin: $(dts[2].time)
        T  Bin: $(dts[3].time)
        R  Bin: $(dts[4].time)
        Z  Bin: $(dts[5].time)
        Y  Bin: $(dts[6].time)
        X  Bin: $(dts[7].time)
    """
end

function getEnergy(v)
    0.5 * me * dot(v, v)
end

function getVelocity(ε)
    sqrt(ε * 2 / me)
end


function adaptiveMeshlog(func, lower, upper, δ)
    # Create an adaptive mesh that has a data point for each "δ" change
    # Start with the lower value
    # Find the func of that value
    # Find the value that has a func value δ different to that of the previous step, this value
    # should be larger
    # stored_xunit = unit(lower)*0+(1*unit(lower))
    stored_yunit = unit(func(lower)) * 0 + (1 * unit(func(lower)))
    # @debug """
    # $(stored_yunit)
    # """
    xs = [lower]
    y = log10(ustrip(func(xs[end])))
    n = 0
    # @debug """
    # $(xs)
    # $y
    # """
    while xs[end] < upper
        n += 1
        n > 100 && break
        sol = find_zeros((z) -> func(z) - 10^(y * (1 + δ)) * stored_yunit, xs[end], upper)
        if isempty(sol)
            sol = find_zeros((z) -> func(z) - 10^(y * (1 - δ)) * stored_yunit, xs[end], upper)
            # @debug """
            # Adaptive Mesh Down
            # $(xs[end])
            # $(y)
            # $(sol)
            # """
            y *= 1 - δ
        else
            # @debug """
            # Adaptive Mesh Up
            # $(xs[end])
            # $(y)
            # $(sol)
            # """
            y *= 1 + δ
        end
        isempty(sol) && break
        push!(xs, minimum(sol))
    end
    push!(xs, upper)

    return xs
end

function adaptiveMesh(func, lower, upper, δ)
    # Create an adaptive mesh that has a data point for each "δ" change
    # Start with the lower value
    # Find the func of that value
    # Find the value that has a func value δ different to that of the previous step, this value
    # should be larger
    xs = [lower]
    y = func(xs[end])
    n = 0
    while xs[end] < upper
        n += 1
        n > 100 && break
        sol = find_zeros((z) -> func(z) - y * (1 + δ), xs[end], upper)
        if isempty(sol)
            sol = find_zeros((z) -> func(z) - y * (1 - δ), xs[end], upper)
            @debug """
            Adaptive Mesh Down
            $(xs[end])
            $(y)
            $(sol)
            """
            y *= 1 - δ
        else
            @debug """
            Adaptive Mesh Up
            $(xs[end])
            $(y)
            $(sol)
            """
            y *= 1 + δ
        end
        isempty(sol) && break
        push!(xs, minimum(sol))
    end
    push!(xs, upper)

    return xs
end

@generated function structprint(this::T)::T where {T}
    expr = Expr(:block)
    contents_expr = Expr(:tuple)
    for fname in fieldnames(this)
        push!(contents_expr.args, :(print(this.$fname)))
    end

    push!(expr.args, :(contents = $(contents_expr)))
    push!(expr.args, :(obj = T(contents...)))

    expr
end

end
