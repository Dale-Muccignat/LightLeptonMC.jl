module MC

using MCModels
using ..Init
using PhysicalConstants.CODATA2018, Unitful, AutoParameters, Random, StaticArrays, LinearAlgebra

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

function flight(particle,model,params,verbose)
    meas = initMeas(params)

    alive = true
    while alive
        verbose && @info "Particle energy before: " getEnergy(particle)
        # if verbose @info "Particle velocity before: " particle.v end
        # if verbose @info "Particle position before: " particle.r end
        verbose && @info "Frequency: " getTotalFreq(particle,model,params,verbose)

        dt = getDt(particle,model,params,verbose)

        updateParticle!(dt.time,particle,model.En0,model.n0,params,verbose)

        event!(dt,particle,model,params,meas,verbose)

        verbose && @info "Particle energy after: " getEnergy(particle)
        # if verbose @info "Particle velocity after: " particle.v end
        # if verbose @info "Particle position after: " particle.r end

        alive = checkAlive(particle,model,params,verbose)
    end
    return meas
end

function checkAlive(particle,model,params,verbose)
    return particle.t < model.t_f
end


#### * Collisions

function event!(dt::dt_COLLISION,particle,model,params,meas,verbose)
    verbose && @info "Energy Bin collision: " dt
    coll_type = getCollType(particle,model,params,verbose)

    if coll_type == :Null
    elseif coll_type == :ELASTIC
        elasticColl!(getRelVel(particle,model),particle,model,params,meas,verbose)
    elseif coll_type == :EXCITATION
    elseif coll_type == :IONIZATION
    elseif coll_type == :ATTACHMENT
    end
end

function getRelVel(particle,model)
    if model.T != 0.
        # u, g = findu(particle,vars,consts,g,u,w,electron_v2eV) #TODO findu
    else
        # g = copy(particle.v)
    end
    return g = copy(particle.v)
end


function getCollType(particle,model,params,verbose)
    R = rand()
    CS = 0.
    en = getEnergy(particle)

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
    # meas.no_elas += 1
    particle.v = v2
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

function getDtEBinsPos(En0::UniformField,particle,model,params)
	v_perp = particle.v[1]^2 + particle.v[2]^2 #Seperate out of the vz velocity
	if (particle.v[3] < 0.) && ((params.ebins[particle.ebin]/electron_v2eV - v_perp) > 0.) && (particle.ebin!= 1)
		#If the velocity is less than zero, i.e if the energy bin needs to go down
		#The second condition is to ensure that the energy can go down considering the perpindicular velocity... or something along those lines. It just makes sure the squareroot doesn't return an imaginary number
		vf_z = sqrt(params.ebins[particle.ebin]/electron_v2eV - v_perp)
		# vf_z is the velocity the particle will need to reach the lower boundry of the energy bin
		return dt_ENERGY(-(particle.v[3]+vf_z) /params.a[3], Down())
		#Just finding the time until the particles velocity reaches vf_z
	else#if (particle.ebin_no < length(vars.ebins))
		# here the vz is positive, and also making sure we aren't at the max bin value
		#Don't need to check the second  condition above as the particle can only increase in energy
		vf_z = sqrt(params.ebins[particle.ebin+1]/electron_v2eV - v_perp)
		#same as before except with the lower boundry of the next bin up
		return dt_ENERGY(-(particle.v[3]-vf_z) /params.a[3], Up())
	end
    return dt_ENERGY(Inf, Up())
end

function event!(dt::dt_ENERGY,particle,model,params,meas,verbose)
    verbose && @info "Energy Bin Change: " dt
    if typeof(dt.direction) == Up
        particle.ebin += 1
    elseif typeof(dt.direction) == Down
        particle.ebin -= 1
    end
end

function updateParticle!(dt,particle,En0::UniformField,n0::UniformDensity,params,verbose)
    particle.r += (particle.v + 0.5 * params.a * dt)*dt
    particle.v += params.a * dt
    particle.t += dt
end

# * Fetch Functions

function getTotalFreq(particle,model,params,verbose)::Float64

    params.total_freq[particle.ebin]

end

function getDensity(particle,n0::UniformDensity,verbose)::Float64
    n0.val
end

function getDt(particle,model,params,verbose)
    dt_coll = dt_COLLISION(randexp()/(getTotalFreq(particle,model,params,verbose)*getDensity(particle,model.n0,verbose)))
    dt_energy = getDtEBinsPos(model.En0,particle,model,params)
    return findMinTime((dt_coll,dt_energy),verbose=verbose)
end

function findMinTime(kw;verbose)

    verbose && @info [k.time for k in kw]
    dt = minimum([k.time for k in kw])
    # kw[sortperm([k.time for k in kw]::Vector{Float64})[1]]
    for k in kw
        if dt == k.time
            return k
        end
    end
    
    return kw[1]

end


end
