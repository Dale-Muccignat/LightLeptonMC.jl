module Utils

using MCSimUnits, LinearAlgebra, Roots
export debugParticle, getEnergy, getVelocity, structprint, adaptiveMesh

function debugParticle(header,particle)
    @debug """
    Particle $header: 
      ε: $(uconvert(eV,getEnergy(particle.v)))
    pos: x=$(particle.r[1]) y=$(particle.r[2]) z=$(particle.r[3]) 
    vel: x=$(particle.v[1]) y=$(particle.v[2]) z=$(particle.v[3]) 
    rad: $(sqrt(particle.r[2]^2+particle.r[3]^2))
    """
end

function getEnergy(v)
    0.5*me*dot(v,v)
end

function getVelocity(ε)
    sqrt(ε*2/me)
end

function adaptiveMesh(func,lower,upper,δ)
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
        sol = find_zeros((z) -> func(z)-y*(1+δ),xs[end],upper)
        if isempty(sol)
            sol = find_zeros((z) -> func(z)-y*(1-δ),xs[end],upper)
            @debug """
            Adaptive Mesh Down
            $(xs[end])
            $(y)
            $(sol)
            """
            y *= 1-δ
        else
            @debug """
            Adaptive Mesh Up
            $(xs[end])
            $(y)
            $(sol)
            """
            y *= 1+δ
        end
        isempty(sol) && break
        push!(xs,minimum(sol))
    end
    push!(xs,upper)
     
    return xs
end

@generated function structprint(this::T)::T where T
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
