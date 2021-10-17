module MC

using ..Structs
using ..Init
using PhysicalConstants.CODATA2018, Unitful, AutoParameters, Random, StaticArrays, LinearAlgebra, AndExport

export flight

const amu = ustrip(float(AtomicMassConstant))
const q = ustrip(float(ElementaryCharge))
const me = ustrip(float(ElectronMass))
const k = ustrip(float(BoltzmannConstant))
const ħ = ustrip(float(ReducedPlanckConstant))

const electron_v2eV = 0.5 * me / q
#TODO: Keep the units in

abstract type Direction end
abstract type Event end
struct Down <: Direction end
struct Up <: Direction end

struct dt_EVENT
    time::Float64 
    type::Symbol 
    direction::Symbol
end


function flight(particle,model,params,verbose)
    meas_inner = initMeas(model,params)

    alive = true
    while alive
        verbose && @info "Particle energy before: " getEnergy(particle.v)/params.energy_scale

        dt = getDt(model.type,particle,model,params,verbose)

        updateParticle!(dt.time,particle,model.En0,model.n0,params,verbose)

        event!(dt,particle,model,params,meas_inner,verbose)

        verbose && @info "Particle energy after: " getEnergy(particle.v)/params.energy_scale

        alive = checkAlive(model.type,particle,model,params,verbose)
        # if verbose && meas_inner.measg.coll_count.elastic > 5 alive = false end
    end
    return meas_inner
end

function checkAlive(type::SWARM,particle,model,params,verbose)
    return particle.t < model.type.t_f*params.time_scale
end

function checkAlive(type::LμJ,particle,model,params,verbose)

    if getEnergy(particle.v) < type.thresh*params.energy_scale
        verbose && @info "Absorbed"
    end

    if particle.rbin == length(params.rbins)
        verbose && @info "Detected"
    end
    particle.rbin != length(params.rbins) && getEnergy(particle.v) > type.thresh*params.energy_scale
end

function event!(dt,particle,model,params,meas_inner,verbose)
    verbose && @info "==================================="
    if dt.type == :dt_collision
        verbose && @info "Collision"
        eventCollision!(dt,particle,model,params,meas_inner,verbose)
    elseif dt.type == :dt_energy
        verbose && @info "Energy"
        eventEnergy!(dt,particle,model,params,meas_inner,verbose)
    elseif dt.type == :dt_tbin
        verbose && @info "TBin"
        eventTbin!(dt,particle,model,params,meas_inner,verbose)
    elseif dt.type == :dt_rbin
        verbose && @info "RBin"
        eventRBin!(dt,particle,model,params,meas_inner,verbose)
    elseif dt.type == :dt_zbin
        verbose && @info "ZBin"
        eventZBin!(dt,particle,model,params,meas_inner,verbose)
    elseif dt.type == :dt_ybin
        verbose && @info "YBin"
        eventYBin!(dt,particle,model,params,meas_inner,verbose)
    else
        
    end
end

@xport @fastmath function updateParticle!(dt,particle,En0::UniformField,n0::UniformDensity,params,verbose)
    particle.r += (particle.v + 0.5 * params.a * params.n0 * dt)*dt
    particle.v += params.a * params.n0 * dt
    particle.t += dt
end

@xport @fastmath function updateParticle!(dt,particle,En0,n0::LμJDensity,params,verbose)
    particle.r += particle.v * dt
    particle.t += dt
end
#### * TBins Measurements
function eventTbin!(dt,particle,model,params,meas_inner,verbose)
    if meas_inner.meast.r != nothing 
        @inbounds @. meas_inner.meast.v²[particle.tbin,:] += particle.v^2 * particle.weight
        @inbounds @. meas_inner.meast.r²[particle.tbin,:] += particle.r^2 * particle.weight
        @inbounds @. meas_inner.meast.r[particle.tbin,:] += particle.r * particle.weight
        @inbounds @. meas_inner.meast.v[particle.tbin,:] += particle.v * particle.weight
        @inbounds @. meas_inner.meast.rv[particle.tbin,:] += particle.v * particle.r * particle.weight
        @inbounds meas_inner.meast.w[particle.tbin] += particle.weight
        @inbounds meas_inner.meast.εᵢ[particle.tbin] += dot(particle.v,particle.v) * particle.weight
    end

    meas_inner.measLμJ.zy_bins != nothing ? meas_inner.measLμJ.zy_bins[particle.zbin,particle.ybin,particle.tbin] += particle.weight : nothing
    
    particle.tbin += 1
end


#### * Collisions

function eventCollision!(dt::dt_EVENT,particle,model,params,meas_inner,verbose)
    verbose && @info "Collision: " dt
    cs_ind = getCollType(particle,model,params,verbose)

    verbose && @info "Collision Type: " params.splat_gas.type[cs_ind]

    if cs_ind == 0
    elseif params.splat_gas.type[cs_ind] == :ELASTIC
        elasticColl!(cs_ind,getRelVel(cs_ind,particle,model,params),particle,model,params,meas_inner,verbose)
    elseif params.splat_gas.type[cs_ind] == :EXCITATION
        excitationColl!(cs_ind,getRelVel(cs_ind,particle,model,params),particle,model,params,meas_inner,verbose)
    elseif params.splat_gas.type[cs_ind] == :IONIZATION
        ionisationColl!(cs_ind,model.ion_type,getRelVel(cs_ind,particle,model,params),particle,model,params,meas_inner,verbose)
    elseif params.splat_gas.type[cs_ind] == :ATTACHMENT
        attachmentColl!(cs_ind,getRelVel(cs_ind,particle,model,params),particle,model,params,meas_inner,verbose)
    end

    updateEbin!(particle,params)
end

function getRelVel(cs_no,particle,model,params)
    if model.T.val == 0.
        g = copy(particle.v)
    else
        g = findu(cs_no,particle,model,params)
    end
    return g 
end


@fastmath function getCollType(particle,model,params,verbose) @inbounds begin
    R = rand()
    CS = 0.
    en = getEnergy(particle.v)/params.energy_scale
    v = sqrt(dot(particle.v,particle.v))/params.length_scale

    for cs_ind in eachindex(params.splat_gas.cs)
        @inbounds CS += params.splat_gas.cs[cs_ind](en)*1e-20 * params.splat_gas.w[cs_ind] * v
        # @show CS
        # @show getTotalFreq(particle,model,params,verbose)

        if R < CS/getTotalFreq(particle,model,params)
            return cs_ind
        end
    end
    return 0
end end

@fastmath function elasticColl!(cs_no,g,particle,model,params,meas_inner,verbose)
    coschi, sinchi, phi = scatteringAngle(g,particle,model,params,verbose)
    v2 = sqrt(dot(g,g)) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
    v2 = shift(v2,g)
    v2 = (me*particle.v + params.splat_gas.m[cs_no]* (g-particle.v + v2)) /(me + params.splat_gas.m[cs_no])
    meas_inner.measg.coschi += dot(particle.v,v2)/(norm(v2)*norm(particle.v))
    meas_inner.measg.coll_count.elastic += 1
    # meas_inner.no_elas += 1
    particle.v = v2
end

@fastmath function excitationColl!(cs_no,g,particle,model,params,meas_inner,verbose)

    if dot(g,g)-2*params.splat_gas.thresh[cs_no]/me > 0.
        coschi, sinchi, phi = scatteringAngle(g,particle,model,params,verbose)
        v2 = sqrt(dot(g,g)-2*params.splat_gas.thresh[cs_no]/me) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
        v2 = shift(v2,g)
        v2 = (me*particle.v + params.splat_gas.m[cs_no]* (g-particle.v + v2)) /(me + params.splat_gas.m[cs_no])
        meas_inner.measg.coschi += dot(particle.v,v2)/(norm(v2)*norm(particle.v))
        meas_inner.measg.coll_count.excitation += 1
        particle.v = v2
    end
    
end

@fastmath function ionisationColl!(cs_no,ion_type,g,particle,model,params,meas_inner,verbose)

    if dot(g,g)-2*params.splat_gas.thresh[cs_no]/me > 0.
        part_energy, other_part_energy = getIonisationEnergy(params.splat_gas.thresh[cs_no],ion_type,g)

        coschi, sinchi, phi = scatteringAngle(g,particle,model,params,verbose)

        v2 = sqrt(2*part_energy/me) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
        v2 = shift(v2,g)
        v2 = (me*particle.v + params.splat_gas.m[cs_no]* (g-particle.v + v2)) /(me + params.splat_gas.m[cs_no])

        particle.weight *= 2.
        meas_inner.measg.coschi += dot(particle.v,v2)/(norm(v2)*norm(particle.v))
        meas_inner.measg.coll_count.ionisation += 1
        particle.v = v2
    end
    
end

function getIonisationEnergy(thresh,ion_type::UniformIonisation,g)
    total_energy = getEnergy(g)*q - thresh
    other_part_energy = total_energy * rand()
    part_energy = total_energy - other_part_energy 
    return part_energy, other_part_energy
end
function getIonisationEnergy(thresh,ion_type::RatioIonisation,g)
    total_energy = getEnergy(g)*q - thresh
    if rand() > 0.5 
        other_part_energy = total_energy * ion_type.ratio
    else
        other_part_energy = total_energy * (1-ion_type.ratio)
    end
    part_energy = total_energy - other_part_energy 
    return part_energy, other_part_energy
end

function attachmentColl!(cs_no,g,particle,model,params,meas_inner,verbose)
    particle.weight /= 2.
end

function scatteringAngle(g,particle,model,params,verbose)
    coschi = 1-2*rand()
    sinchi = sqrt(1-coschi^2)
    phi = rand()*2*pi
    return coschi, sinchi, phi
end

@fastmath function shift(v_new::SArray{Tuple{3},Float64,1,3},v::SArray{Tuple{3},Float64,1,3})::SArray{Tuple{3},Float64,1,3}
    relvel_mag = norm(v)
    orig_costheta = v[3] / relvel_mag
    orig_sintheta = sqrt(1 - orig_costheta^2)
    if orig_sintheta == 0
        orig_cosphi = 1.
        orig_sinphi = 0.
    else
        orig_cosphi = v[1] / relvel_mag / orig_sintheta
        orig_sinphi = v[2] / relvel_mag / orig_sintheta
    end
    
    Ry = @SMatrix [orig_costheta 0. orig_sintheta ; 0. 1. 0. ; -orig_sintheta 0. orig_costheta]
    Rz = @SMatrix [orig_cosphi -orig_sinphi 0. ; orig_sinphi orig_cosphi 0. ; 0. 0. 1.]

	return Rz * Ry * v_new
end
#### * Update Bin numbers

function updateEbin!(particle,params)
    ea = getEnergy(particle.v)
	while ea < params.ebins[particle.ebin]
		particle.ebin-= 1
	end
    if particle.ebin != length(params.ebins)
        while ea > params.ebins[particle.ebin+1]
            particle.ebin+= 1
        end
    end
    
end

function eventEnergy!(dt::dt_EVENT,particle,model,params,meas_inner,verbose)
    verbose && @info "Energy Bin Change: " dt
    if dt.direction == :Up
        particle.ebin += 1
    elseif dt.direction == :Down
        particle.ebin -= 1
    end
end

function eventRBin!(dt::dt_EVENT,particle,model,params,meas_inner,verbose)
    verbose && @info "RBin Bin Change: " dt
    if dt.direction == :Up
        particle.rbin += 1
    elseif dt.direction == :Down
        particle.rbin -= 1
    end

    if particle.rbin == length(params.rbins)
        meas_inner.measLμJ.EELS[findMinInd(params.ebinsLμJ,getEnergy(particle.v))] += particle.weight
    end
    
end

function eventYBin!(dt::dt_EVENT,particle,model,params,meas_inner,verbose)
    verbose && @info "yBin Bin Change: " dt
    if dt.direction == :Up
        particle.ybin += 1
    elseif dt.direction == :Down
        particle.ybin -= 1
    end
end

function eventZBin!(dt::dt_EVENT,particle,model,params,meas_inner,verbose)
    verbose && @info "zBin Bin Change: " dt
    if dt.direction == :Up
        particle.zbin += 1
    elseif dt.direction == :Down
        particle.zbin -= 1
    end
end

# * Get Bin Times

function getDtEBinsPos(En0::UniformField,n0::UniformDensity,particle,model,params)
	v_perp = particle.v[1]^2 + particle.v[2]^2 #Seperate out of the vz velocity
	if (particle.v[3] < 0.) && ((params.ebins[particle.ebin]/electron_v2eV - v_perp) > 0.) && (particle.ebin!= 1)
		#If the velocity is less than zero, i.e if the energy bin needs to go down
		#The second condition is to ensure that the energy can go down considering the perpindicular velocity... or something along those lines. It just makes sure the squareroot doesn't return an imaginary number
		vf_z = sqrt(params.ebins[particle.ebin]/electron_v2eV - v_perp)
		# vf_z is the velocity the particle will need to reach the lower boundry of the energy bin
		return dt_EVENT(-(particle.v[3]+vf_z) /(params.a[3] * params.n0),:dt_energy, :Down)
		#Just finding the time until the particles velocity reaches vf_z
	else#if (particle.ebin_no < length(vars.ebins))
		# here the vz is positive, and also making sure we aren't at the max bin value
		#Don't need to check the second  condition above as the particle can only increase in energy
		vf_z = sqrt(params.ebins[particle.ebin+1]/electron_v2eV - v_perp)
		#same as before except with the lower boundry of the next bin up
		return dt_EVENT(-(particle.v[3]-vf_z) /(params.a[3] * params.n0), :dt_energy,:Up)
	end
    return dt_ENERGY(Inf, :dt_energy, Up())
end

@fastmath function getRDt(particle,model,params)
	#NOTE: Only for fields of 0
    p = particle.r
    v = particle.v
    r21 = params.rbins[particle.rbin]
    r22 = params.rbins[particle.rbin+1]
	up, down = Inf, Inf
	a = v[3]^2+v[2]^2
	b = 2*(p[3]*v[3]+v[2]*p[2])
	c = p[3]^2+p[2]^2-r21+1e-6
    delta = b^2-4*a*c
	if (delta > 0.) #&& (c > 0.)
		t1 = (-b+sqrt(delta))/(2*a)
		t2 = (-b-sqrt(delta))/(2*a)
		if t1 <= 0.
			t1 = Inf
		end
		if t2 <= 0.
			t2 = Inf
		end
		down = min(t1,t2)
		if r21 == 0.
			down = Inf
		end
	end
	c = c+r21-r22-2e-6
    delta = b^2-4*a*c
	if (delta > 0.) #&& (c < 0.)
		t1 = (-b+sqrt(delta))/(2*a)
		t2 = (-b-sqrt(delta))/(2*a)
		if t1 <= 0.
			t1 = Inf
		end
		if t2 <= 0.
			t2 = Inf
		end
		up = min(t1,t2)
	end
    if up == min(up,down)
        return dt_EVENT(up,:dt_rbin,:Up)
    else
        return dt_EVENT(down,:dt_rbin,:Down)
    end
end

function getCDt(particle,model,params)
    dt_coll = dt_EVENT(randexp()/(getTotalFreq(particle,model,params)*getDensity(particle,model.n0,params)*params.length_scale^3)::Float64,:dt_collision,:Up)
end

function getZDt(particle,model,params)
    if params.zbins != nothing
        if (particle.v[3] < 0.) && (particle.zbin != 1)
            return dt_EVENT((params.zbins[particle.zbin]-1e-12-particle.r[3])/particle.v[3], :dt_zbin,:Down)
        elseif (particle.zbin < length(params.zbins)) && (particle.v[3] > 0.)
            return dt_EVENT((params.zbins[particle.zbin+1]+1e-12-particle.r[3])/particle.v[3], :dt_zbin, :Up)
        end
    end
    return dt_EVENT(Inf,:dt_zbin,:Up)
end

function getYDt(particle,model,params)
    if params.ybins != nothing
        if (particle.v[2] < 0.) && (particle.ybin != 1)
            return dt_EVENT((params.ybins[particle.ybin]-1e-12-particle.r[2])/particle.v[2], :dt_ybin,:Down)
        elseif (particle.ybin < length(params.ybins)) && (particle.v[2] > 0.)
            return dt_EVENT((params.ybins[particle.ybin+1]+1e-12-particle.r[2])/particle.v[2], :dt_ybin, :Up)
        end
    end
    return dt_EVENT(Inf,:dt_ybin,:Up)
end

function getTDt(particle,model,params)
    if params.tbins != nothing 
        if particle.tbin < length(params.tbins) && particle.t > params.tbins[1]
            return dt_tbin = dt_EVENT(params.tbins[particle.tbin+1] - particle.t, :dt_tbin,:Up)
        end
    end
    
    return dt_EVENT(Inf,:dt_tbin,:Up)

end

# * Fetch Functions

@xport function getTotalFreq(particle,model,params)
    params.total_freq[particle.ebin]
end

@xport function getDensity(particle,n0::UniformDensity,params)
    params.n0
end
@xport function getDensity(particle,n0::LμJDensity,params)
    params.n0[particle.rbin]
end

@xport function getDt(type::SWARM,particle,model,params,verbose)
    verbose && @info "Current time: " particle.t
    dt_coll = getCDt(particle,model,params)
    dt_energy = getDtEBinsPos(model.En0,model.n0,particle,model,params)
    dt_tbin = getTDt(particle,model,params)
    dts = (dt_coll,dt_energy,dt_tbin)
    verbose && @info dts
    dt_times = (dt_coll.time,dt_energy.time,dt_tbin.time)
    return dts[findMinTime(dt_times,verbose=verbose)]
end



@xport function getDt(type::LμJ,particle,model,params,verbose)
    verbose && @info "Current time: " particle.t
    # @info particle.r ./ params.length_scale
    dt_coll = getCDt(particle,model,params)
    dt_rbin = getRDt(particle,model,params)
    dt_tbin = getTDt(particle,model,params)
    dt_zbin = getZDt(particle,model,params)
    dt_ybin = getYDt(particle,model,params)
    dts = (dt_coll,dt_rbin,dt_zbin,dt_ybin,dt_tbin)
    verbose && @info dts
    dt_times = (dt_coll.time,dt_rbin.time,dt_zbin.time,dt_ybin.time,dt_tbin.time)
    return dts[findMinTime(dt_times,verbose=verbose)]
end

function findMinTime(times;verbose)

    verbose && @info times
    dt = minimum(times)
    for i in 1:length(times)
        if dt == times[i]
            return i
        end
    end
    
    return 1

end


function findTotalCF(en,model,params)
    CS = 0.
    for cs_ind in eachindex(params.splat_gas.cs)
        @inbounds CS += params.splat_gas.cs[cs_ind](en) * params.splat_gas.w[cs_ind] * sqrt(en/electron_v2eV)
    end

    return CS
end
#### * Temperature


@fastmath function findu(cs_ind,particle,model,params)
    #TODO: make generic function based on inputs
	w = sqrt(2*k*model.T.val/params.splat_gas.m[cs_ind])
    v = particle.v*params.time_scale/params.length_scale
    # @show getEnergy(v)

	he = electron_v2eV * (norm(v) + 5. *w)^2 
	le = max( 1e-25 , electron_v2eV * (norm(v) - 5. *w).^2 )
    # @show he
    # @show le

    #TODO: get actual CF from CS
    
	# @inbounds hcf = params.splat_gas.cs[cs_ind](he) * sqrt(he/electron_v2eV)
	# @inbounds lcf = params.splat_gas.cs[cs_ind](le) * sqrt(le/electron_v2eV)
	@inbounds hcf = findTotalCF(he,model,params)
	@inbounds lcf = findTotalCF(le,model,params)


    g = copy(v)
	mcf = max(hcf,lcf)
	while true
		# usually it only goes through one while loop iteration
		u = SArray{Tuple{3},Float64,1,3}(randn(3))*w./sqrt(2.) # could only run this
		g = v .- u
        # @show v
        # @show u
        # @show g
		eg = electron_v2eV * dot(g,g) 
        # @show eg
		@inbounds ccf = findTotalCF(eg,model,params)
        # @show lcf
        # @show hcf
        # @show ccf
        # @show mcf

        if ccf/mcf == 0.
            throw("Zero relative energy")
        end
        
		if rand() < ccf/mcf
			break
		end
	end
	return g*params.length_scale/params.time_scale
end


end
