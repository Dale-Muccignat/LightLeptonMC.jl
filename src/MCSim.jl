"""
MCSim

This package takes gases from Gases.jl and performs a Monte Carlo simulation.

`MC(::Gas,::Params)`

"""
module MCSim


using Reexport, Distributed
@reexport using Gases
@reexport using MCModels
export mc, MCINIT, MCPARAMS

include("Init.jl")
@reexport using .Init

include("MC.jl")
@reexport using .MC

include("Measurements.jl")
@reexport using .Measurements

function mc(model=MCMODEL(),params=initParams(MCMODEL());verbose=false)

    meas = initMeas(params)

    meas = @distributed (+) for n in 1:model.n
        if n%(model.n*0.01) == 0.
            @info "Progress: "*string(100*n/model.n)*"%"
        end

        flight(PART(model,params),model,params,verbose)
    end
    @show meas
    return process(meas,params,model);
    

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
Base.:+(x::MEAST,y::MEAST) = GenericAddition(x,y)::MEAST
# Base.:+(x::MEAS,y::MEAS) = GenericAddition(x,y)::MEAS
Base.:+(x::MEASS,y::MEASS) = GenericAddition(x,y)::MEASS
Base.:+(x::MEASC,y::MEASC) = GenericAddition(x,y)::MEASC
Base.:+(x::COLLCOUNT,y::COLLCOUNT) = GenericAddition(x,y)::COLLCOUNT
function Base.:+(x::MEAS,y::MEAS)
    return MEAS(GenericAddition(x.measg,y.measg),GenericAddition(x.meast,y.meast),GenericAddition(x.meass,y.meass),GenericAddition(x.measc,y.measc))
end


end
