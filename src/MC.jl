module MC

using MCModels
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

@AutoParm struct dt_COLLISION <: Event
    time::Float64 = Inf
end

@AutoParm struct dt_ENERGY <: Event
    time::Float64 = Inf
    direction::Direction = Up()
end

@AutoParm mutable struct dt_EVENT
    time::Float64 = Int
    type::Symbol = :dt_collision
    direction::Symbol = :Up
end


function flight(particle,model,params,verbose)
    meas = initMeas(params)

    alive = true
    while alive
        verbose && @info "Particle energy before: " getEnergy(particle.v)
        # if verbose @info "Particle velocity before: " particle.v end
        # if verbose @info "Particle position before: " particle.r end
        verbose && @info "Frequency: " getTotalFreq(particle,model,params,verbose)

        dt = getDt(particle,model,params,verbose)

        updateParticle!(dt.time,particle,model.En0,model.n0,params,verbose)

        event!(dt,particle,model,params,meas,verbose)

        verbose && @info "Particle energy after: " getEnergy(particle.v)
        # if verbose @info "Particle velocity after: " particle.v end
        # if verbose @info "Particle position after: " particle.r end

        alive = checkAlive(particle,model,params,verbose)
    end
    return meas
end

function checkAlive(particle,model,params,verbose)
    return particle.t < model.t_f
end

function event!(dt,particle,model,params,meas,verbose)
    if dt.type == :dt_collision
        eventCollision!(dt,particle,model,params,meas,verbose)
    elseif dt.type == :dt_energy
        eventEnergy!(dt,particle,model,params,meas,verbose)
    elseif dt.type == :dt_tbin
        eventTbin!(dt,particle,model,params,meas,verbose)
    else
        
    end
end

#### * TBins Measurements
function eventTbin!(dt,particle,model,params,meas,verbose)

    @inbounds @. meas.meast.v²[particle.tbin,:] += particle.v^2 * particle.weight
    @inbounds @. meas.meast.r²[particle.tbin,:] += particle.r^2 * particle.weight
    @inbounds @. meas.meast.r[particle.tbin,:] += particle.r * particle.weight
    @inbounds @. meas.meast.v[particle.tbin,:] += particle.v * particle.weight
    @inbounds @. meas.meast.rv[particle.tbin,:] += particle.v * particle.r * particle.weight
    @inbounds meas.meast.w[particle.tbin] += particle.weight
    @inbounds meas.meast.εᵢ[particle.tbin] += dot(particle.v,particle.v) * particle.weight
    particle.t += dt.time
    particle.tbin += 1
end


#### * Collisions

function eventCollision!(dt::dt_EVENT,particle,model,params,meas,verbose)
    verbose && @info "Energy Bin collision: " dt
    coll_type = getCollType(particle,model,params,verbose)

    if coll_type == :Null
    elseif coll_type == :ELASTIC
        elasticColl!(getRelVel(particle,model),particle,model,params,meas,verbose)
    elseif coll_type == :EXCITATION
    elseif coll_type == :IONIZATION
    elseif coll_type == :ATTACHMENT
    end

    updateEbin!(particle,params)
end

function getRelVel(particle,model)
    if model.T != 0.
        # u, g = findu(particle,vars,consts,g,u,w,electron_v2eV) #TODO findu
    else
        # g = copy(particle.v)
    end
    return g = copy(particle.v)
end


function getCollType(particle,model,params,verbose)::Symbol
    R = rand()
    CS = 0.
    en = getEnergy(particle.v)

    for (cs,type) in splatCS(model)
        CS += cs(en)*1e-20 * getDensity(particle,model.n0,verbose) * sqrt(2*q*en/me)
        if R < CS/getTotalFreq(particle,model,params,verbose)
            return type
        end
    end
    return :Null
end

function elasticColl!(g,particle,model,params,meas,verbose)
    verbose && @info "Elastic Collision"
    coschi, sinchi, phi = scatteringAngle(g,particle,model,params,verbose)
    v2 = sqrt(dot(g,g)) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
    v2 = shift(v2,g)
    v2 = (me*particle.v + model.gas.Mₑₗ[1]*amu * (g-particle.v + v2)) /(me + model.gas.Mₑₗ[1]*amu)
    meas.measg.coschi += dot(particle.v,v2)/(norm(v2)*norm(particle.v))
    meas.measg.coll_count.elastic += 1
    # meas.no_elas += 1
    particle.v = v2
end

function excitationColl!(g,particle,model,params,meas,verbose)
    verbose && @info "Excitation Collision"
    #TODO: Add Thresholds to cross sections

    if dot(g,g)-2*model.gas.thresh/me > 0.
        v2 = sqrt(dot(g,g)-2*vars.thresh/me) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
        v2 = shift(v2,g)
        v2 = (me*particle.v + model.gas.Mₑₗ[1]*amu * (g-particle.v + v2)) /(me + model.gas.Mₑₗ[1]*amu)
        meas.measg.coschi += dot(particle.v,v2)/(norm(v2)*norm(particle.v))
        meas.measg.coll_count.inelastic += 1
        particle.v = v2
    end
    
end

# function ionisationColl!(iontype::UniformIonisation,g,particle,model,params,meas,verbose)
#     verbose && @info "Elastic Collision"

#     total_energy = getEnergy(particle.v)*q - vars.thresh
#     other_part_energy = total_energy * rand()

#     part_energy = total_energy - other_part_energy 

#     coschi, sinchi, phi = scatteringAngle(vars,type_no,particle.ebin_no)

#     v2 = sqrt(2*part_energy/consts.me) * SArray{Tuple{3},Float64,1,3}(cos(phi)*sinchi,sin(phi)*sinchi,coschi)
#     v2 = shift(v2,g)
#     v2 = (me*particle.v + model.gas.Mₑₗ[1]*amu * (g-particle.v + v2)) /(me + model.gas.Mₑₗ[1]*amu)

#     meas.measg.coschi += dot(particle.v,v2)/(norm(v2)*norm(particle.v))
#     meas.measg.coll_count.ionisation += 1
#     particle.v = v2
# end

function updateEbin!(particle,params)
    ea = getEnergy(particle.v)
	while ea < params.ebins[particle.ebin]
		# @show "no"
		particle.ebin-= 1
	end
	while ea > params.ebins[particle.ebin+1]
		# @show "no"
		particle.ebin+= 1
	end
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
# * Measurements

# * Bin Events

function getDtEBinsPos(En0::UniformField,n0::UniformDensity,particle,model,params)
	v_perp = particle.v[1]^2 + particle.v[2]^2 #Seperate out of the vz velocity
	if (particle.v[3] < 0.) && ((params.ebins[particle.ebin]/electron_v2eV - v_perp) > 0.) && (particle.ebin!= 1)
		#If the velocity is less than zero, i.e if the energy bin needs to go down
		#The second condition is to ensure that the energy can go down considering the perpindicular velocity... or something along those lines. It just makes sure the squareroot doesn't return an imaginary number
		vf_z = sqrt(params.ebins[particle.ebin]/electron_v2eV - v_perp)
		# vf_z is the velocity the particle will need to reach the lower boundry of the energy bin
		return dt_EVENT(-(particle.v[3]+vf_z) /(params.a[3] * n0.val),:dt_energy, :Down)
		#Just finding the time until the particles velocity reaches vf_z
	else#if (particle.ebin_no < length(vars.ebins))
		# here the vz is positive, and also making sure we aren't at the max bin value
		#Don't need to check the second  condition above as the particle can only increase in energy
		vf_z = sqrt(params.ebins[particle.ebin+1]/electron_v2eV - v_perp)
		#same as before except with the lower boundry of the next bin up
		return dt_EVENT(-(particle.v[3]-vf_z) /(params.a[3] * n0.val), :dt_energy,:Up)
	end
    return dt_ENERGY(Inf, Up())
end

function eventEnergy!(dt::dt_EVENT,particle,model,params,meas,verbose)
    verbose && @info "Energy Bin Change: " dt
    if dt.direction == :Up
        particle.ebin += 1
    elseif dt.direction == :Down
        particle.ebin -= 1
    end
end

@xport function updateParticle!(dt,particle,En0::UniformField,n0::UniformDensity,params,verbose)
    particle.r += (particle.v + 0.5 * params.a * n0.val * dt)*dt
    particle.v += params.a * n0.val * dt
    particle.t += dt
end

# * Fetch Functions

function getTotalFreq(particle,model,params,verbose)::Float64
    params.total_freq[particle.ebin]
end

function getDensity(particle,n0::UniformDensity,verbose)::Float64
    n0.val
end

@xport function getDt(particle,model,params,verbose)
    verbose && @info "Current time: " particle.t
    dt_coll = dt_EVENT(randexp()/(getTotalFreq(particle,model,params,verbose)*getDensity(particle,model.n0,verbose)),:dt_collision,:Up)
    dt_energy = getDtEBinsPos(model.En0,model.n0,particle,model,params)
    dt_tbin = dt_EVENT(params.tbins[particle.tbin+1] - particle.t, :dt_tbin,:Up)
    dts = (dt_coll,dt_energy,dt_tbin)
    dt_times = (dt_coll.time,dt_energy.time,dt_tbin.time)
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


end
