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

    meas = initMeas(model,params)

    # @show params.rbins
    # meas = @distributed (+) for n in 1:model.n
    for n in 1:model.n

        if n%(model.n*0.1) == 0.
            @info "Progress: "*string(100*n/model.n)*"%"
        end

        meas += flight(PART(model.type,model,params),model,params,verbose)
    end
    # @show params.rbins
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
