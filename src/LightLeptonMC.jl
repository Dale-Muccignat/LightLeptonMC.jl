"""
LightLeptonMC

This package takes gases from LightLeptonTools.jl and performs a Monte Carlo simulation in either a swarm or LμJ configuration. An example swarm simulation is included below.

```julia
meas, swarm, params = mc(MaxwellModel();E=1.,T=0.,N=1000,S=0.,verbose=false)
```

For more information, type ?mc
"""
module LightLeptonMC
#TODO: Add liquid
#TODO: Add Potential
#TODO: Add Density Gradient (this is probably best solved through Danny's code though)
#TODO: Reform struct system into a generic Detector/Object system

using LightLeptonTools
using Distributed, Statistics, Unitful, Optim, AutoParameters, StaticArrays, LinearAlgebra, Random, Measurements, SparseArrays, QuadGK, Interpolations
using DataFramesMeta, DataFrames

Unitful.register(@__MODULE__)
function __init__()
    Unitful.register(@__MODULE__)
    return Unitful.preferunits(m̃, s̃)
end

const electron_v2eV = 0.5 * me

include("Structs.jl")
export MCMODEL, initParams, CELL, LμJ, InitVelocity, InitSpatial, UniformDensity, UniformField, GaussianTemporal, UniformTemp

include("Utils.jl")

include("Init.jl")

include("TimeOfFlight.jl")

include("Collisions.jl")

include("Events.jl")

include("MC.jl")

include("Processing.jl")

export mc, MCPARAMS, saveBatch, saveDist

macro jld2save(filename, vars...)
    # Copy paste
    # writeexprs = Vector{Expr}(undef, length(vars))
    # for i = 1:length(vars)
    #     writeexprs[i] = :(JLD2.write(f, $(string(vars[i])), $(esc(vars[i])), wsession))
    # end

    JLD2 = esc(:JLD2)
    writeexprs = map(vars) do var
        :($JLD2.write(f, $(string(var)), $(esc(var)), wsession))
    end

    return quote
        if !isdefined(@__MODULE__, :JLD2)
            error("Need to have the JLD2 module included for this to work.")
        end

        $JLD2.jldopen($(esc(filename)), true, true, true, Base.IOStream) do f
            wsession = $JLD2.JLDWriteSession()
            $(Expr(:block, writeexprs...))
        end
    end
end

@doc "This function takes the output of prepareSwarmData and duplicates it for a given gas"
function mc(gas::LLCS, prepared_swarm_data::DataFrame; kwargs...)
    mirrored_swarm_data = @transform prepared_swarm_data @byrow begin
        :ys = mirrorSwarmDataRow(gas, :x, :T, :P, :types; kwargs...)
        :gas = gas
    end
    return mirrored_swarm_data
end

@doc "For a given  row, this function takes the output of prepareSwarmData and duplicates it for a given gas"
function mirrorSwarmDataRow(
        gas::LLCS,
        Es::Vector{<:ReducedVoltage},
        T::Temperature,
        P::Union{Missing, Pressure},
        types::Vector{Symbol};
        kwargs...
    )
    #TODO: Better Convergence handling...
    sols = [
        mc(gas; P, T, E, kwargs...)[2]
            for E in Es
    ]
    return [[t == :ε ? s.ε : t == :W ? s.bW : t == :WF ? s.fW : t == :nDLF ? s.fDl : t == :nDTF ? s.fDt : t == :nDL ? s.bDl : t == :Keff ? s.R : t == :nDT ? s.bDt : error("$(t) not implemented") for s in sols] for t in types]
    #TODO: rename these so that it's a 1:1 conversion to FVB
end

function mc(gas::LLCS; P::Pressure = 0kPa, S::Energy = 0.01eV, E::ReducedVoltage = 1.0Td, T::Temperature = 0.0K, N = 1000, t::Time = Inf * s, meas = [:t], ion_type = UniformIonisation(:EXPLICIT), n0::NumberDensity = 1.0e23 / m^3)
    if t == Inf * s
        t = findSimTime(gas, S, E, T, n0)[1] / n0
    end
    model = MCMODEL(meas = meas, type = SWARM(t), n = N, En0 = UniformField(E), T = UniformTemp(T), ion_type = ion_type, v_0 = InitVelocity(S, [0, 0, 1]), n0 = UniformDensity(n0), optimise_spatial = false)
    return mc(model, initParams(gas, model))
end

function findSimTime(gas, S, E, T, n0)
    #= curr_time = 1e-9s * n0 =#
    curr_time = 1.0e13s / m^3
    unconverged = true
    N::Int64 = 10000
    N_og = N
    iterations = 0
    measures_o = nothing
    swarm_o = nothing
    while unconverged
        if iterations > 100 #|| N > 1000000
            throw("Convergence Failed. iters = $(iterations), N = $(N)")
        end
        iterations += 1
        @debug "Current Time: " curr_time
        model = MCMODEL(type = SWARM(curr_time / n0), n = N, n0 = UniformDensity(n0), En0 = UniformField(E), T = UniformTemp(T), ion_type = UniformIonisation(:EXPLICIT), v_0 = InitVelocity(S, [0, 0, 1]))
        params = initParams(gas, model)
        measures, swarm = mc(model, params)
        particle_change = mean(measures.meast.w[(3 * end ÷ 4):end])
        first, second = mean(measures.meast.εᵢ[(end ÷ 3):(2 * end ÷ 3)]), mean(measures.meast.εᵢ[(2 * end ÷ 3):end])
        first_std, second_std = std(measures.meast.εᵢ[(end ÷ 3):(2 * end ÷ 3)]), std(measures.meast.εᵢ[(2 * end ÷ 3):end])
        if particle_change < N_og * 0.75
            # if we've lost half our particles
            N = round(Int, N * N_og / particle_change)
            @debug "Lost $(N_og - particle_change) particles, changing N to $(N)"
            #= iterations += 1 =#
            continue # retry
        elseif N > 100 && particle_change > N_og * 2.0
            N = round(Int, N * N_og / particle_change)
            if N < 100
                N = 100
            end
            @debug "Gained $(particle_change - N_og) particles, changing N to $(N)"
            #= iterations += 1 =#
            continue # retry
        end
        #= @info "STDs: " first_std second_std third_std =#
        @debug "Averages: $(first) $(second)"
        @debug "Diff: $(abs(second - first) / second)"
        @debug "Particles: $(particle_change)"
        if abs(second - first) / second > 0.01
            curr_time *= 1.25
        else
            #= curr_time *= 2.0 # For good measure =#
            @debug "Chose: " curr_time
            unconverged = false
            measures_o = measures
            swarm_o = swarm
        end
    end
    return curr_time, N_og / N, measures_o, swarm_o
end


function mc(model::MCMODEL, gas::LLCS)
    return mc(model, initParams(gas, model))
end

function mc(model::MCMODEL, gases, indicies, w = 1.0, ref_gas = MaxwellModel(), save = false, dcs = false)
    distribution_data::Dict{Int64, Tuple{Vector{Float64}, Vector{Float64}}} = Dict()
    # @info "Batch: " indicies[1]
    for i in indicies
        # @info "Number: " i
        if w == 1.0
            params = initParams(gases[i], model)
        else
            params = initParams(gases[i] * w + ref_gas * (1 .- w), model)
        end

        meas, swarm, params = mc(model, params)
        if dcs
            distribution_data[i] = (params.χbins[4:(end - 1)], sum(meas.measLμJ.εχ[:, :], dims = 1)'[4:(end - 1)])
        else
            distribution_data[i] = (params.ebins, meas.measLμJ.ε)
        end
    end
    if save
        saveBatch(distribution_data, model, indicies, w)
    end

    return distribution_data
end


"""

General models
```julia
meas, swarm, params = mc(model::MCMODEL,params::MCPARAMS)
mean, swarm, params = mc(model,gas)
```

Swarm models
```julia
meas, swarm, params = mc(MaxwellModel();E=1.,T=0.,N=1000,S=0.)
```

"""
function mc(model::MCMODEL, params::MCPARAMS)
    # @debug model.n=1

    #scale!(model,model.type,params)
    #' meas_outer = initMeas(model, params)
    #' meass = [deepcopy(meas_outer) for _ in 1:Threads.nthreads()]
    @debug "Measurements initialised."
    #scale!(params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    N = model.n
    if N / nprocs() < 1000
        batches = [n for n in Iterators.partition(1:N, ceil(Int64, N / nprocs()))]
    else
        batches = [n for n in Iterators.partition(1:N, 100)]
    end

    meas = @distributed (+) for n in eachindex(batches)
        Ni = length(batches)
        #' Threads.@threads for n in 1:N
        #' id = Threads.threadid()
        #' @show id
        #' for n in eachindex(batches)
        #' meas = @distributed (+) for n in 1:model.n
        #=     @debug "Distributing..." =#
        #==#
        if n % (Ni * 0.1) == 0.0
            @info "Progress: " * string(100 * n / Ni) * "%"
        end

        # TODO: Account for no of ionisation events
        # meas_inner = flight(PART(model,params),model,params,verbose)
        # if model.ion_type.method == :EXPLICIT

        # end
        # [meas_inner]
        #= [flight(PART(model, params), model, params)] =#
        #' flight!([PART(model, params)], model, params, meass[id])
        flight([PART(model, params) for _ in batches[n]], model, params)
    end
    @debug "Finished Sim"
    # @show params.rbins
    #unscale!(params)
    #unscale!(model.type,model,params,meas[1])
    #' meas = sum(meass)
    return process(model.type, meas, params, model)
end

function mcGetParticles(model::MCMODEL, params::MCPARAMS, coll_limit = Inf)
    # @debug model.n=1

    #scale!(model,model.type,params)
    #= meas = initMeas(model, params) =#
    @debug "Measurements initialised."
    #scale!(params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    particles = @distributed (vcat) for n in 1:model.n
        @debug "Distributing..."

        if n % (model.n * 0.1) == 0.0
            @info "Progress: " * string(100 * n / model.n) * "%"
        end

        # TODO: Account for no of ionisation events
        # meas_inner = flight(PART(model,params),model,params,verbose)
        # if model.ion_type.method == :EXPLICIT

        # end
        # [meas_inner]
        flightGetParticles(PART(model, params), model, params, coll_limit = coll_limit)
    end
    @debug "Finished Sim"
    # @show params.rbins
    #unscale!(params)
    return particles
end

function mcInitParticles(particles_og::Vector{<:PART}, model::MCMODEL, params::MCPARAMS; N = 0, meas = initMeas(model, params))
    N != 0 && (model.n += N)
    N == 0 && (N = model.n)
    if N / nprocs() < 1000
        batches = [n for n in Iterators.partition(1:N, ceil(Int64, N / nprocs() / 4))]
    else
        batches = [n for n in Iterators.partition(1:N, 100)]
    end
    #= batches = [n for n in Iterators.partition(1:N, 1)] =#

    meas_new = @distributed (+) for n in eachindex(batches)
        @debug "Distributing..."

        #= if n % (model.n * 0.1) == 0.0 =#
        #=     @info "Progress: " * string(100 * n / model.n) * "%" =#
        #= end =#

        # TODO: Account for no of ionisation events
        # meas_inner = flight(PART(model,params),model,params,verbose)
        # if model.ion_type.method == :EXPLICIT

        # end
        # [meas_inner]
        particles = [PART(model, params) for _ in batches[n]]
        sample_particles = sample(particles_og, length(batches[n]))
        #= particle.r = sample_particle.r =#
        for i in eachindex(batches[n])
            particles[i].v = sample_particles[i].v
            particles[i].weight = sample_particles[i].weight
            particles[i].cbin = findMinInd(params.cbins, getEnergy(particles[i].v))
            if params.ebins != nothing
                particles[i].ebin = findMinInd(params.ebins, getEnergy(particles[i].v))
            end
            if params.zbins != nothing
                particles[i].zbin = findMinInd(params.zbins, particles[i].r[3])
            end
        end
        [flight(particles, model, params)]
    end
    @debug "Finished Sim"
    # @show params.rbins
    #unscale!(params)
    #unscale!(model.type,model,params,meas[1])
    meas = unweights(meas, model)
    meas += meas_new[1]
    return process(model.type, meas, params, model)
end

# The goal of the function is to provide a framework to progressibely build the electron swarm. There are two facets. I need to ensure that there are a sufficient number of particles and that they have to simulated for long enough.

function mcAddParticles(model, params; N = 0, meas = initMeas(model, params))
    # this takes a vector of particles and simulates them until some end condition is met...
    # Status: work in progress. The flaw is having to redefine params everytime anyway.
    N != 0 && (model.n += N)
    N == 0 && (N = model.n)
    @debug N
    if N / nprocs() < 1000
        batches = [n for n in Iterators.partition(1:N, ceil(Int64, N / nprocs() / 4))]
    else
        batches = [n for n in Iterators.partition(1:N, 100)]
    end
    #= batches = [n for n in Iterators.partition(1:N,100)] =#
    #= batches = [[n] for n in 1:N] =#

    meas_new = @distributed (+) for n in eachindex(batches)
        @debug "Distributing..."

        if n % (N * 0.1) == 0.0
            @info "Progress: " * string(100 * n / N) * "%"
        end

        # TODO: Account for no of ionisation events
        # meas_inner = flight(PART(model,params),model,params,verbose)
        # if model.ion_type.method == :EXPLICIT

        # end
        # [meas_inner]
        [flight([PART(model, params) for _ in batches[n]], model, params)]
    end
    #= return meas_new =#
    meas = unweights(meas, model)
    meas += meas_new[1]
    @debug "Finished Sim"
    return process(model.type, meas, params, model)
end

function buildModel()
    println("Welcome to the model builder helper")
    println("What type of model [LuJ/Swarm/Cell]:")
    modeltype = readline()
    while modeltype ∉ ["LuJ", "Swarm", "Cell"]
        println("Error: Model not found")
        println("What type of model [LuJ/Swarm/Cell]:")
        modeltype = readline()
    end

    println("How many particles?")
    ns = readline()
    n = 0
    while n == 0
        try
            n = parse(Int64, ns)
        catch
            println("Error: Not an Int64")
            println("How many particles?")
            ns = readline()
        end
    end

    println("What initial energy (in eV)?")
    εs = readline()
    ε = Inf
    while ε == Inf
        try
            ε = parse(Float64, εs)
        catch
            println("Error: Not an Float64")
            println("What initial energy (in eV)?")
            εs = readline()
        end
    end
    return
end

function GenericAddition!(this::T, other::T) where {T}
    for f in fieldnames(T)
        a = getfield(this, f)
        b = getfield(other, f)
        if isa(a, AbstractArray) && isa(b, AbstractArray)
            a .+= b
        else
            setfield!(this, f, a + b)
        end
    end
    return this
end

#' @generated function GenericAddition(this::T, other::T)::T where {T}
#'     expr = Expr(:block)
#'     contents_expr = Expr(:tuple)
#'     for fname in fieldnames(this)
#'         push!(contents_expr.args, :(this.$fname + other.$fname))
#'     end
#'
#'     push!(expr.args, :(contents = $(contents_expr)))
#'     push!(expr.args, :(obj = T(contents...)))
#'
#'     expr
#' end

Base.:+(x::Nothing, y::Nothing) = nothing
Base.:*(x::Nothing, y::Float64) = nothing
Base.:/(x::Nothing, y::Float64) = nothing
Base.:*(x::Nothing, y) = nothing
Base.:/(x::Nothing, y) = nothing
Unitful.upreferred(x::Nothing) = nothing
Unitful.uconvert(y, x::Nothing) = nothing

Base.:+(x::SWARMMEAS, y::SWARMMEAS) = GenericAddition!(x, y)::SWARMMEAS
Base.:+(x::MEAST, y::MEAST) = GenericAddition!(x, y)::MEAST
Base.:+(x::MEASLμJ, y::MEASLμJ) = GenericAddition!(x, y)::MEASLμJ
Base.:+(x::MEASCELL, y::MEASCELL) = GenericAddition!(x, y)::MEASCELL
Base.:+(x::MEASC, y::MEASC) = GenericAddition!(x, y)::MEASC
Base.:+(x::MEASS, y::MEASS) = GenericAddition!(x, y)::MEASS
function Base.:+(x::MEAS, y::MEAS)
    return MEAS(GenericAddition!(x.measLμJ, y.measLμJ), GenericAddition!(x.measCell, y.measCell), GenericAddition!(x.measg, y.measg), GenericAddition!(x.meast, y.meast), GenericAddition!(x.meass, y.meass), GenericAddition!(x.measc, y.measc))
end

function saveBatch(dists, model, indicies, w)
    i = indicies[1]
    j = indicies[end]
    en = ustrip(model.v_0.ε)
    n = model.n
    filename = "ML/" * "$w" * "/" * "$en" * "eV/" * "$en" * "_no" * "$n" * "_#" * "$i" * ":" * "$j" * "_" * string(Dates.today()) * "_" * string(Dates.hour(now())) * "-" * string(Dates.minute(now())) * "-" * string(Dates.second(now())) * "-" * string(Dates.millisecond(now()))
    @jld2save filename * "TEMP" * ".jld2" dists
    cp(filename * "TEMP" * ".jld2", filename * ".jld2")
    return rm(filename * "TEMP" * ".jld2")
end

function saveDist(name, model, params, meas)
    en = ustrip(model.v_0.ε)
    n = model.n
    ebins = params.ebins
    dist = meas.measLμJ.ε
    filename = "$name" * "_" * "$en" * "_no" * "$n" * "_" * string(Dates.today()) * "_" * string(Dates.hour(now())) * "-" * string(Dates.minute(now())) * "-" * string(Dates.second(now())) * "-" * string(Dates.millisecond(now()))
    @jld2save filename * "TEMP" * ".jld2" dist ebins
    cp(filename * "TEMP" * ".jld2", filename * ".jld2")
    return rm(filename * "TEMP" * ".jld2")
end

end
