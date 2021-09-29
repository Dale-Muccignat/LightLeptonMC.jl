"""
MCLuJ

This package takes gases from Gases.jl and performs a Monte Carlo simulation.

`MC(::Gas,::Params)`

"""
module MCLuJ

using Reexport
include("Init.jl")
@reexport using .Init

end
