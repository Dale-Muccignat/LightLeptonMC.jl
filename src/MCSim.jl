"""
MCSim

This package takes gases from Gases.jl and performs a Monte Carlo simulation in either a swarm or LμJ configuration. An example swarm simulation is included below.

```julia
    meas, swarm = mc(MaxwellModel();E=1.,T=0.,N=1000,S=0.,verbose=false)
    mean, swarm = mc(model,gas)
```

"""
module MCSim
#TODO: Add Attachment, super-eleastic, liquid
#TODO: Add Potential
#TODO: Add Density Gradient (this is probably best solved through Danny's code though)

using Reexport, Distributed, JLD2, Dates, AndExport
@reexport using Gases
@reexport using FVB
export mc, MCINIT, MCPARAMS, saveBatch, saveDist

include("Structs.jl")
@reexport using .Structs

include("Init.jl")
@reexport using .Init

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

function mc(gas::Gas;S=0.,E=1.,T=0.,N=1000,t=1e-5,verbose=false)
    model = MCMODEL(type=SWARM(t),n=N,En0=UniformField(E),T=UniformTemp(T),v_0=InitVelocity(S,[0,0,1]))
    mc(model,initParams(gas,model),verbose=verbose)
end

function mc(model::MCMODEL,gas::Gas;verbose=false)
    mc(model,initParams(gas,model),verbose=verbose)
end

function mc(model::MCMODEL,gases::Vector{Gas},indicies;verbose=false,w=1.0,ref_gas=MaxwellModel(),save=false)
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
        distribution_data[i] = (params.ebins,meas.measLμJ.ε)
    end
    if save
        saveBatch(distribution_data,model,indicies,w)
    end
    
    return distribution_data
end


function mc(model::MCMODEL,params::MCPARAMS;verbose=false)

    meas = initMeas(model,params)
    scale!(params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    meas = @distributed (+) for n in 1:model.n

        if n%(model.n*0.1) == 0.
            @info "Progress: "*string(100*n/model.n)*"%"
        end

        meas_inner = flight(PART(model,params),model,params,verbose)
        [meas_inner]
    end
    # @show params.rbins
    unscale!(params)
    unscale!(model.type,params,meas[1])
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

Base.:+(x::SWARMMEAS,y::SWARMMEAS) = GenericAddition(x,y)::SWARMMEAS
Base.:+(x::MEAST,y::MEAST) = GenericAddition(x,y)::MEAST
Base.:+(x::MEASLμJ,y::MEASLμJ) = GenericAddition(x,y)::MEASLμJ
Base.:+(x::MEASCELL,y::MEASCELL) = GenericAddition(x,y)::MEASCELL
Base.:+(x::MEASC,y::MEASC) = GenericAddition(x,y)::MEASC
Base.:+(x::MEASS,y::MEASS) = GenericAddition(x,y)::MEASS
Base.:+(x::COLLCOUNT,y::COLLCOUNT) = GenericAddition(x,y)::COLLCOUNT
function Base.:+(x::MEAS,y::MEAS)
    return MEAS(GenericAddition(x.measLμJ,y.measLμJ),GenericAddition(x.measCell,y.measCell),GenericAddition(x.measg,y.measg),GenericAddition(x.meast,y.meast),GenericAddition(x.meass,y.meass),GenericAddition(x.measc,y.measc))
end

function saveBatch(dists,model,indicies,w)
    i = indicies[1]
    j = indicies[end]
    en = model.v_0.ε
    n = model.n
    filename = "ML/"*"$w"*"/"*"$en"*"eV/"*"$en"*"_no"*"$n"*"_#"*"$i"*":"*"$j"*"_" * string(Dates.today()) * "_" * string(Dates.hour(now())) * "-" * string(Dates.minute(now())) * "-" * string(Dates.second(now())) * "-" * string(Dates.millisecond(now()))
    @jld2save filename*"TEMP"*".jld2" dists
    cp(filename*"TEMP"*".jld2",filename*".jld2")
    rm(filename*"TEMP"*".jld2")
end

function saveDist(name,model,params,meas)
    en = model.v_0.ε
    n = model.n
    ebins = params.ebins
    dist = meas.measLμJ.ε
    filename = "$name"*"_"*"$en"*"_no"*"$n"*"_" * string(Dates.today()) * "_" * string(Dates.hour(now())) * "-" * string(Dates.minute(now())) * "-" * string(Dates.second(now())) * "-" * string(Dates.millisecond(now()))
    @jld2save filename*"TEMP"*".jld2" dist ebins
    cp(filename*"TEMP"*".jld2",filename*".jld2")
    rm(filename*"TEMP"*".jld2")
end

end
