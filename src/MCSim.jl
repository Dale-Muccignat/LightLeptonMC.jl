"""
MCSim

This package takes gases from Gases.jl and performs a Monte Carlo simulation in either a swarm or LμJ configuration. An example swarm simulation is included below.

```julia
    meas, swarm, params = mc(MaxwellModel();E=1.,T=0.,N=1000,S=0.,verbose=false)
```

For more information, type ?mc
"""
module MCSim
#TODO: Add super-eleastic, liquid
#TODO: Add Potential
#TODO: Add Density Gradient (this is probably best solved through Danny's code though)

using Reexport, Distributed, JLD2, Dates, AndExport, Statistics, Unitful
using Gases
using FVB
export mc, MCINIT, MCPARAMS, saveBatch, saveDist
@reexport using MCSimUnits


include("Structs.jl")
@reexport using .Structs

include("Utils.jl")
@reexport using .Utils

include("Init.jl")
@reexport using .Init

include("TimeOfFlight.jl")
@reexport using .TimeOfFlight

include("MC.jl")
@reexport using .MC

include("Measurements.jl")
@reexport using .Measurements


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

    quote
        if !isdefined(@__MODULE__, :JLD2)
            error("Need to have the JLD2 module included for this to work.")
        end

        $JLD2.jldopen($(esc(filename)), true, true, true, Base.IOStream) do f
            wsession = $JLD2.JLDWriteSession()
            $(Expr(:block, writeexprs...))
        end
    end
end



function mc(gas::Gas;S=1e-5,E=1.,T=0.,N=1000,t=Inf,meas=[:t],ion_type=RatioIonisation(),n0=1e23,verbose=false)
    if t == Inf
        t = findSimTime(gas,S,E,T,meas,verbose)
    end

    model = MCMODEL(meas=meas,type=SWARM(t*s),n=N,En0=UniformField(E*Td),T=UniformTemp(T*K),ion_type=ion_type,v_0=InitVelocity(S*eV,[0,0,1]),n0=UniformDensity(n0/m^3))
    mc(model,initParams(gas,model),verbose=verbose)
end

function findSimTime(gas,S,E,T,meas,verbose)
    curr_time = 1e-7
    unconverged = true
    N = 100
    iterations = 0
    while unconverged
        iterations += 1
        verbose && @info "Current Time: " curr_time
        model = MCMODEL(meas=meas,type=SWARM(curr_time*s),n=1000,En0=UniformField(E*Td),T=UniformTemp(T*K),ion_type=UniformIonisation(),v_0=InitVelocity(S*eV,[0,0,1]))
        measures, swarm, params = mc(model,initParams(gas,model))
        first, second = mean(measures.meast.εᵢ[1:end÷2]),mean(measures.meast.εᵢ[end÷2:end])
        verbose && @info "Averages: " first second
        verbose && @info "Diff: " abs(second-first)/second
        if abs(second-first)/second > 0.01
            curr_time *= 2
        else
            curr_time *= 2 # For good measure
            @info "Chose: " curr_time
            unconverged = false
        end
        if iterations > 10
            @info "Not Converged"
            throw("Convergence Failed.")
        end
    end
    return curr_time
end

function mc(model::MCMODEL,gas::Gas;verbose=false)
    mc(model,initParams(gas,model),verbose=verbose)
end

function mc(model::MCMODEL,gases,indicies;verbose=false,w=1.0,ref_gas=MaxwellModel(),save=false,dcs=false)
    distribution_data::Dict{Int64,Tuple{Vector{Float64},Vector{Float64}}} = Dict()
    # @info "Batch: " indicies[1]
    for i in indicies
        # @info "Number: " i
        if w == 1.
            params = initParams(gases[i],model)
        else
            params = initParams(gases[i]*w + ref_gas*(1 .-w),model)
        end

        meas, swarm, params= mc(model,params,verbose=verbose)
        if dcs
            distribution_data[i] = (params.χbins[4:end-1],sum(meas.measLμJ.εχ[:,:],dims=1)'[4:end-1])
        else
            distribution_data[i] = (params.ebins,meas.measLμJ.ε)
        end
    end
    if save
        saveBatch(distribution_data,model,indicies,w)
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
    meas, swarm, params = mc(MaxwellModel();E=1.,T=0.,N=1000,S=0.,verbose=false)
```

"""
function mc(model::MCMODEL,params::MCPARAMS;verbose=false)
    # @debug model.n=1

    #scale!(model,model.type,params)
    meas = initMeas(model,params)
    @debug "Measurements initialised."
    #scale!(params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    meas = @distributed (+) for n in 1:model.n
        @debug "Distributing..."

        if n%(model.n*0.1) == 0.
            @info "Progress: "*string(100*n/model.n)*"%"
        end

        # TODO: Account for no of ionisation events
        # meas_inner = flight(PART(model,params),model,params,verbose)
        # if model.ion_type.method == :EXPLICIT

        # end
        # [meas_inner]
        [flight(PART(model,params),model,params,verbose)]
    end
    @debug "Finished Sim"
    # @show params.rbins
    #unscale!(params)
    #unscale!(model.type,model,params,meas[1])
    return process(model.type,meas[1],params,model)..., params;

end

@xport function buildModel()
    println("Welcome to the model builder helper")
    println("What type of model [LuJ/Swarm/Cell]:")
    modeltype = readline()
    while modeltype ∉ ["LuJ","Swarm","Cell"]
        println("Error: Model not found")
        println("What type of model [LuJ/Swarm/Cell]:")
        modeltype = readline()
    end

    println("How many particles?")
    ns = readline()
    n=0
    while n == 0
        try
            n=parse(Int64,ns)
        catch
            println("Error: Not an Int64")
            println("How many particles?")
            ns = readline()
        end
    end

    println("What initial energy (in eV)?")
    εs = readline()
    ε=Inf
    while ε == Inf
        try
            ε=parse(Float64,εs)
        catch
            println("Error: Not an Float64")
            println("What initial energy (in eV)?")
            εs = readline()
        end
    end
end


@generated function GenericAddition(this::T,other::T)::T where T
    expr = Expr(:block)
    contents_expr = Expr(:tuple)
    for fname in fieldnames(this)
        push!(contents_expr.args, :(this.$fname + other.$fname))
    end

    push!(expr.args, :(contents = $(contents_expr)))
    push!(expr.args, :(obj = T(contents...)))

    expr
end

Base.:+(x::Nothing,y::Nothing) = nothing
Base.:*(x::Nothing,y::Float64) = nothing
Base.:/(x::Nothing,y::Float64) = nothing
Base.:*(x::Nothing,y) = nothing
Base.:/(x::Nothing,y) = nothing
Unitful.upreferred(x::Nothing) = nothing
Unitful.uconvert(y,x::Nothing) = nothing

Base.:+(x::SWARMMEAS,y::SWARMMEAS) = GenericAddition(x,y)::SWARMMEAS
Base.:+(x::MEAST,y::MEAST) = GenericAddition(x,y)::MEAST
Base.:+(x::MEASLμJ,y::MEASLμJ) = GenericAddition(x,y)::MEASLμJ
Base.:+(x::MEASCELL,y::MEASCELL) = GenericAddition(x,y)::MEASCELL
Base.:+(x::MEASC,y::MEASC) = GenericAddition(x,y)::MEASC
Base.:+(x::MEASS,y::MEASS) = GenericAddition(x,y)::MEASS
function Base.:+(x::MEAS,y::MEAS)
    return MEAS(GenericAddition(x.measLμJ,y.measLμJ),GenericAddition(x.measCell,y.measCell),GenericAddition(x.measg,y.measg),GenericAddition(x.meast,y.meast),GenericAddition(x.meass,y.meass),GenericAddition(x.measc,y.measc))
end

function saveBatch(dists,model,indicies,w)
    i = indicies[1]
    j = indicies[end]
    en = ustrip(model.v_0.ε)
    n = model.n
    filename = "ML/"*"$w"*"/"*"$en"*"eV/"*"$en"*"_no"*"$n"*"_#"*"$i"*":"*"$j"*"_" * string(Dates.today()) * "_" * string(Dates.hour(now())) * "-" * string(Dates.minute(now())) * "-" * string(Dates.second(now())) * "-" * string(Dates.millisecond(now()))
    @jld2save filename*"TEMP"*".jld2" dists
    cp(filename*"TEMP"*".jld2",filename*".jld2")
    rm(filename*"TEMP"*".jld2")
end

function saveDist(name,model,params,meas)
    en = ustrip(model.v_0.ε)
    n = model.n
    ebins = params.ebins
    dist = meas.measLμJ.ε
    filename = "$name"*"_"*"$en"*"_no"*"$n"*"_" * string(Dates.today()) * "_" * string(Dates.hour(now())) * "-" * string(Dates.minute(now())) * "-" * string(Dates.second(now())) * "-" * string(Dates.millisecond(now()))
    @jld2save filename*"TEMP"*".jld2" dist ebins
    cp(filename*"TEMP"*".jld2",filename*".jld2")
    rm(filename*"TEMP"*".jld2")
end

end
