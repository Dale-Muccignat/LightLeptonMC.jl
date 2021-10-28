module Structs

using AutoParameters, Gases, StaticArrays, AndExport

export MCMODEL, MCINIT, MCPARAMS, PART, MEAS, MEAST, MEASS, MEASC, findMinInd

abstract type Velocity end
abstract type Spacial end
abstract type Density end
abstract type Field end
abstract type Temp end

@xport struct ZeroSpacial <: Spacial end
@xport @AutoParm struct InitSpacial <: Spacial
    val::AUTO
    axis::AUTO
end
@xport @AutoParm struct GaussianLμJSpacial <: Velocity
    val::AUTO = -7.6e-6
    Δy_source::AUTO = 7.5e-6
    gaussian_axis::AUTO = [0,1,0]
    axis::AUTO = [0,0,1]
end

@xport @AutoParm struct UniformDensity <: Density
    val::AUTO = 1e23
end
@xport @AutoParm struct LμJDensity <: Density
    val::AUTO = 1e28
end
@xport @AutoParm struct ZeroField <: Field
    val::AUTO = 0.
end
@xport @AutoParm struct UniformField <: Field
    val::AUTO = 1.
end
@xport @AutoParm struct UniformTemp <: Temp
    val::AUTO = 0.
end

abstract type Ionisation end

@xport @AutoParm mutable struct UniformIonisation <: Ionisation end
@xport @AutoParm mutable struct RatioIonisation <: Ionisation
    ratio::AUTO = 0.5
end

@xport struct ZeroVelocity <: Velocity end

@xport @AutoParm struct InitVelocity <: Velocity
    ε::AUTO = 100.
    axis::AUTO = [0,0,1]
end
@xport @AutoParm struct GaussianVelocity <: Velocity
    ε::AUTO = 10.
    Δε_source::AUTO = 0.5
    axis::AUTO = [0,0,1]
end

@AutoParm mutable struct MODEL
    E::AUTO = UniformField()
    T::AUTO = 0.
    init::AUTO = MCINIT()
end

@AutoParm mutable struct MCMODEL
    type::AUTO = SWARM()
    n::Int64 = 100
    v_0::AUTO = ZeroVelocity()
    r_0::AUTO = ZeroSpacial()
    n0::AUTO = UniformDensity()
    En0::AUTO = UniformField()
    T::AUTO = UniformTemp()
    ion_type::AUTO = UniformIonisation()
    meas::AUTO = [:t]
end

@xport @AutoParm mutable struct SWARM
    t_f::AUTO = 2e-5
end

@xport @AutoParm mutable struct LμJ
    Δε_detector::AUTO = 0.9
    thresh::AUTO = 1.
    jet_width::AUTO = 15e-6
end

@xport function LμJ(;Δε_detector=0.9,thresh=1.,jet_width=15e-6)
    LμJ(Δε_detector,thresh,jet_width)
end



@AutoParm mutable struct MCPARAMS
    cbins::AUTO = [0.,(10. .^(-3:0.01:2))...]
    ebins::AUTO = nothing
    gas::AUTO = nothing
    splat_gas::AUTO = nothing
    total_freq::AUTO = nothing
    a::AUTO = nothing
    # n0::Union{Vector{Float64},Float64,Nothing} = nothing
    n0::AUTO = nothing
    rbins::Union{Vector{Float64},Nothing} = nothing
    χbins::Union{Vector{Float64},Nothing} = nothing
    tbins::Union{Vector{Float64},Nothing} = collect(LinRange(0,1e-5,10))
    zbins::Union{Vector{Float64},Nothing} = nothing
    xbins::Union{Vector{Float64},Nothing} = nothing
    ybins::Union{Vector{Float64},Nothing} = nothing
    length_scale::AUTO = nothing
    time_scale::AUTO = nothing
    energy_scale::AUTO = nothing

end

function MCPARAMS(;cbins=nothing,ebins=nothing,gas=nothing,splat_gas=nothing,total_freq=nothing,a=nothing,n0=nothing,rbins=nothing,χbins=nothing,tbins=nothing,zbins=nothing,xbins=nothing,ybins=nothing,length_scale=noting,time_scale=nothing,energy_scale=nothing)
    MCPARAMS(cbins,ebins,gas,splat_gas,total_freq,a,n0,rbins,χbins,tbins,zbins,xbins,ybins,length_scale,time_scale,energy_scale)
end


@xport @AutoParm mutable struct PART
	v::SArray{Tuple{3},Float64,1,3} = [0.,0.,0.]
	r::SArray{Tuple{3},Float64,1,3} = [0.,0.,0.]
	t::Float64 = 0.
	weight::Float64 = 1.
	tbin::Int64 = 1
	cbin::Int64 = 1
	zbin::Int64 = 1
	ybin::Int64 = 1
    xbin::Int64 = 1
    rbin::Int64 = 1
end

@AutoParm mutable struct MEAS
    measLμJ::AUTO = MEASLμJ()
    measg::AUTO = MEASG()
    meast::AUTO = MEAST()
    meass::AUTO = MEASS()
    measc::AUTO = MEASC()
end

function MEAS(type::SWARM,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(),MEASG(),MEAST(params),MEASS(),MEASC())
end
function MEAS(type::LμJ,model,params::MCPARAMS=MCPARAMS())
    MEAS(MEASLμJ(model,params),MEASG(),MEAST(),MEASS(),MEASC())
end

@xport @AutoParm mutable struct COLLCOUNT
    elastic::AUTO = 0
    excitation::AUTO = 0
    ionisation::AUTO = 0
    attachment::AUTO = 0
end

@xport @AutoParm mutable struct MEASLμJ
    EELS::AUTO = nothing
    χ::AUTO = nothing
    zy_bins::AUTO = nothing
end

function MEASLμJ(model::MCMODEL,params::MCPARAMS)
    EELS,χ,zyt_bins = nothing, nothing, nothing
    if params.ybins != nothing && params.zbins != nothing
        zyt_bins = zeros(length(params.zbins),length(params.ybins),length(params.tbins))
    end
    
    if params.ebins != nothing
        EELS = zeros(length(params.ebins))
    end

    if params.χbins != nothing
        χ = zeros(length(params.χbins))
    end
    
    return MEASLμJ(EELS, χ, zyt_bins)

end

@AutoParm mutable struct MEASG
    coschi::AUTO = 0.
    coll_count::AUTO = COLLCOUNT()
end


@AutoParm mutable struct MEAST
    εᵢ::AUTO = nothing
    r::AUTO = nothing
    v::AUTO = nothing
    r²::AUTO = nothing
    v²::AUTO = nothing
    rv::AUTO = nothing
    w::AUTO = nothing
    zeff::AUTO = nothing
end

function MEAST(params)

    MEAST(zeros(length(params.tbins)-1),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1,3),
          zeros(length(params.tbins)-1),
          zeros(length(params.tbins)-1)
          )

end

@AutoParm mutable struct MEASS

end

@AutoParm mutable struct MEASC

end

@AutoParm @xport mutable struct SWARMMEAS
    ε::Float64 = 0.
    εs::Vector{Float64} = [0.]
    bW::Float64 = 0.
    bWs::Vector{Float64} = [0.]
    bDt::Float64 = 0.
    bDts::Vector{Float64} = [0.]
    bDl::Float64 = 0.
    bDls::Vector{Float64} = [0.]
end

@fastmath function findMinInd(x::Array{Float64,1},y::Float64)::Int64
	left::Int64 = 1
	right::Int64 = length(x)
	d::Float64  = Inf
	dlast::Float64  = -Inf
	c::Int64 = 0
	while d != dlast
		# c = Int(round((left+right)/2))
		c = div(left+right,2)
		dlast = d
		@inbounds d = x[c]
		if d < y
			left = c
		elseif d >= y
			right = c
		end

	end
	if x[right] <= y
		return right
	else
		return left
	end
end

@fastmath function findMinInd(x::Array{Int64,1},y::Int64)::Int64
	left::Int64 = 1
	right::Int64 = length(x)
	d::Float64  = Inf
	dlast::Float64  = -Inf
	c::Int64 = 0
	while d != dlast
		# c = Int(round((left+right)/2))
		c = div(left+right,2)
		dlast = d
		@inbounds d = x[c]
		if d < y
			left = c
		elseif d >= y
			right = c
		end

	end
	if x[right] <= y
		return right
	else
		return left
	end
end
end # module
