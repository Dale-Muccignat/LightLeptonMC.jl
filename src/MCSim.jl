"""
MCSim

This package takes gases from Gases.jl and performs a Monte Carlo simulation in either a swarm or LμJ configuration. An example workflow is included below.

```
    `mc(::Gas;E,T,verbose)`

```

"""
module MCSim
#TODO: Add Attachment, super-eleastic, liquid
#TODO: Add Potential
#TODO: Add Density Gradient (this is probably best solved through Danny's code though)

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

function mc(gas::Gas;S=0.,E=1.,T=0.,N=1000,verbose=false)
    model = MCMODEL(n=N,En0=UniformField(1.),T=UniformTemp(0.),v_0=InitVelocity(S,[0,0,1]))
    mc(model,initParams(model.type,gas,model),verbose=verbose)
end


function mc(model=MCMODEL(),params=initParams(MCMODEL());verbose=false)

    # meas = initMeas(model,params)
    scale!(params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    meas = @distributed (+) for n in 1:model.n

        if n%(model.n*0.1) == 0.
            @info "Progress: "*string(100*n/model.n)*"%"
        end

        flight(PART(model.type,model,params),model,params,verbose)
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
Base.:+(x::SWARMMEAS,y::SWARMMEAS) = GenericAddition(x,y)::SWARMMEAS
Base.:+(x::MEAST,y::MEAST) = GenericAddition(x,y)::MEAST
Base.:+(x::MEASLμJ,y::MEASLμJ) = GenericAddition(x,y)::MEASLμJ
Base.:+(x::MEASC,y::MEASC) = GenericAddition(x,y)::MEASC
Base.:+(x::MEASS,y::MEASS) = GenericAddition(x,y)::MEASS
Base.:+(x::COLLCOUNT,y::COLLCOUNT) = GenericAddition(x,y)::COLLCOUNT
function Base.:+(x::MEAS,y::MEAS)
    return MEAS(GenericAddition(x.measLμJ,y.measLμJ),GenericAddition(x.measg,y.measg),GenericAddition(x.meast,y.meast),GenericAddition(x.meass,y.meass),GenericAddition(x.measc,y.measc))
end


end
