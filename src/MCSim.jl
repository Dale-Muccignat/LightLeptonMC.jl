"""
MCSim

This package takes gases from Gases.jl and performs a Monte Carlo simulation in either a swarm or LμJ configuration. An example workflow is included below.

First, initilise the chosen model:

```
    MCMODEL()                                      # creates a default swarm configuration
    MCMODEL(type={SWARM(),LμJ()},v_0,r_0,n0,En0)   # creates a custom configuration
```

Then, initilise the params object with the chose Gas:

```
    MCPARAMS()                                     # creates parameters with a Maxwell gas
    MCPARAMS({SWARM,LμJ},::Gas,::MCMODEL)          # creates custom parameters with a gas from Gases.jl
```

Finally, run the Monte Carlo simulation:

```
    mc(::MCMODEL,::MCPARAMS)                       # runs a Monte Carlo simluation using the provided model and configuration
```

"""
module MCSim


using Reexport, Distributed
@reexport using Gases
export mc, MCINIT, MCPARAMS

include("Structs.jl")
@reexport using .Structs

include("Init.jl")
@reexport using .Init

include("MC.jl")
@reexport using .MC

include("Measurements.jl")
@reexport using .Measurements

function mc(model=MCMODEL(),params=initParams(MCMODEL());verbose=false)

    meas = initMeas(model,params)
    scale!(params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    for n in 1:model.n

        if n%(model.n*0.1) == 0.
            @info "Progress: "*string(100*n/model.n)*"%"
        end

        meas += flight(PART(model.type,model,params),model,params,verbose)
    end
    # @show params.rbins
    unscale!(params)
    unscale!(model.type,params,meas)
    return process(model.type,meas,params,model);
    

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
Base.:+(x::MEAST,y::MEAST) = GenericAddition(x,y)::MEAST
Base.:+(x::MEASLμJ,y::MEASLμJ) = GenericAddition(x,y)::MEASLμJ
Base.:+(x::MEASC,y::MEASC) = GenericAddition(x,y)::MEASC
Base.:+(x::MEASS,y::MEASS) = GenericAddition(x,y)::MEASS
Base.:+(x::COLLCOUNT,y::COLLCOUNT) = GenericAddition(x,y)::COLLCOUNT
function Base.:+(x::MEAS,y::MEAS)
    return MEAS(GenericAddition(x.measLμJ,y.measLμJ),GenericAddition(x.measg,y.measg),GenericAddition(x.meast,y.meast),GenericAddition(x.meass,y.meass),GenericAddition(x.measc,y.measc))
end


end
