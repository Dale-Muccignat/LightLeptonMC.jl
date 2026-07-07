#### * Collisions

function eventCollision!(dt::dt_EVENT, particle, model, params, meas_inner, new_particles)
    new_particle = nothing
    cs_ind = getCollType(particle, model, params)
    energy_loss = 0.0
    :SP in model.meas && meas_inner.measg.total_coll == 0 ? energy_loss = getEnergy(particle.v) : nothing
    #' @debug """
    #' Collision:
    #' $(cs_ind == 0 ? "NULL" : params.gas[cs_ind][1].type) $(dt.time)
    #' """
    if cs_ind != 0
        cs, weight, mass = params.gas[cs_ind]
        g = getRelVel(cs_ind, cs, mass, particle, model, params)

        if cs.type === :ELASTIC || cs.type === :EFFECTIVE
            elasticColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
        elseif cs.type === :EXCITATION || cs.type === :ROTATION || cs.type === :VIBRATION || cs.type === :DISSOCIATION
            excitationColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
            # energy_loss = cs.energy_loss
        elseif cs.type === :IONIZATION
            ionisationColl!(cs_ind, cs, mass, model.ion_type, g, particle, model, params, meas_inner, new_particles)
            # energy_loss = cs.energy_loss
        elseif cs.type === :ATTACHMENT || cs.type === :THERMAL || cs.type === :FORMATION || cs.type === :ANNIHILATION
            attachmentColl!(particle, model.ion_type.method)
        elseif cs.type === :ELASTIC_MOMENTUM
            elasticMomentumColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
        elseif cs.type === :ELASTIC_ENERGY
            elasticEnergyColl!(cs_ind, cs, mass, g, particle, model, params, meas_inner)
        end
        meas_inner.measg.coll_count[cs_ind] += 1
        meas_inner.measg.total_coll += 1
    end

    :SP in model.meas && meas_inner.measg.total_coll == 1 ? energy_loss -= getEnergy(particle.v) : nothing

    measureCollision!(energy_loss, particle, model, params, meas_inner)
    updateCbin!(particle, params)
    return new_particle
end

function getRelVel(cs_no, cs, mass, particle, model, params)
    if model.T.val == 0.0 * K
        g = copy(particle.v)
    else
        g = findu(cs_no, cs, mass, particle, model, params)
    end
    #' @debug g
    return g
end

function getCollType(particle, model, params)
    begin
        R = rand()
        CS = 0.0Å^2 * sqrt(eV)
        en = uconvert(eV, getEnergy(particle.v)) #/params.energy_scale
        #' v = sqrt(dot(particle.v, particle.v))#/params.length_scale
        total_freq = getTotalFreq(particle, model, params)

        for (cs_ind, (cs, weight, mass)) in enumerate(params.gas)
            if model.ion_type.method === :WEIGHTING && (cs.type === :ATTACHMENT || cs.type === :THERMAL)
                CS += 2 * CollisionFrequency(cs)(sqrt(en)) * weight
            else
                CS += CollisionFrequency(cs)(sqrt(en)) * weight
            end

            if R < CS / total_freq
                return cs_ind
            end
        end
        return 0
    end
end

function elasticColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)
    #me2 = 4amu
    # @debug uconvert(amu,me), uconvert(amu,mass)
    coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
    v2 = norm(g) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
    v2 = shift(v2, g)
    v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)
    meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
    # meas_inner.no_elas += 1
    return particle.v = v2
end

function elasticMomentumColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)
    #me2 = 4amu
    # @debug uconvert(amu,me), uconvert(amu,mass)
    coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
    v2 = norm(g) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
    v2 = shift(v2, g)
    v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)
    meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
    v2 = v2 * norm(particle.v) / norm(v2)
    # meas_inner.no_elas += 1
    return particle.v = v2
end

function elasticEnergyColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)
    #me2 = 4amu
    # @debug uconvert(amu,me), uconvert(amu,mass)
    coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
    v2 = norm(g) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
    v2 = shift(v2, g)
    v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)
    meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
    # meas_inner.no_elas += 1
    return particle.v = particle.v * norm(v2) / norm(particle.v)
end

function excitationColl!(cs_no, cs, mass, g, particle, model, params, meas_inner)

    return if dot(g, g) - 2 * cs.energy_loss / me > 0m^2 / s^2
        coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
        v2 = sqrt(dot(g, g) - 2 * cs.energy_loss / me) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
        v2 = shift(v2, g)
        v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)
        meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
        particle.v = v2
    end

end

function ionisationColl!(cs_no, cs, mass, ion_type, g, particle, model, params, meas_inner, new_particles)

    return if dot(g, g) - 2 * cs.energy_loss / me > 0m^2 / s^2
        part_energy, other_part_energy = getIonisationEnergy(cs.energy_loss, ion_type, g, particle, model, params)

        coschi, sinchi, phi = scatteringAngle(cs_no, cs, mass, g, particle, model, params)

        v2 = sqrt(2 * part_energy / me) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
        v2 = shift(v2, g)
        v2 = (me * particle.v + mass * (g - particle.v + v2)) / (me + mass)

        meas_inner.measg.coschi += dot(particle.v, v2) / (norm(v2) * norm(particle.v))
        particle.v = v2


        if ion_type.method == :WEIGHTING
            particle.type isa Electron ? particle.weight *= 2.0 : nothing
        elseif ion_type.method == :EXPLICIT
            coschi2, sinchi2, phi2 = scatteringAngle(cs_no, cs, mass, g, particle, model, params)
            v2 = sqrt(2 * other_part_energy / me) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi2) * sinchi2, sin(phi2) * sinchi2, coschi2)
            new_particle = deepcopy(particle)
            meas_inner.measg.coschi += dot(new_particle.v, v2) / (norm(v2) * norm(new_particle.v))
            new_particle.v = v2
            # @info debugParticle("other particle",new_particle)
            updateCbin!(new_particle, params)
            push!(new_particles, new_particle)
        else
            error("Unknown non-conservative method in ion_type")
        end
    else
        @error begin
            """
            Ion threshold greater than energy"
            Cross section #: $(cs_no)
            Threshold: $(cs.energy_loss)
            Dot: $(dot(g, g))
            minus: $(2 * cs.energy_loss / me)
            $(debugParticle("Ionisation collision", particle))

            """
        end
    end

end

function getIonisationEnergy(thresh, ion_type::DiscreteIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type isa Positron
        other_part_energy = ion_type.loss #*params.energy_scale*q
    else
        other_part_energy = total_energy - ion_type.loss #*params.energy_scale*q
    end
    part_energy = total_energy - other_part_energy
    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::UniformIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    other_part_energy = total_energy * rand()
    part_energy = total_energy - other_part_energy
    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::RatioIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type isa Positron
        other_part_energy = total_energy * ion_type.ratio(getEnergy(g))::Float64
        # other_part_energy = 0
    else
        other_part_energy = total_energy * (1 - ion_type.ratio(getEnergy(g))::Float64)
    end
    part_energy = total_energy - other_part_energy
    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::FunctionalIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type isa Positron
        other_part_energy = ion_type.func[particle.cbin](rand()) #*params.energy_scale*q
    else
        other_part_energy = total_energy - ion_type.func[particle.cbin](rand()) #*params.energy_scale*q
    end
    part_energy = total_energy - other_part_energy

    if part_energy < 0
        @info "part" part_energy #/params.energy_scale/q
        @info "other" other_part_energy #/params.energy_scale/q
    end

    return part_energy, other_part_energy
end

function getIonisationEnergy(thresh, ion_type::DiscreteFunctionalIonisation, g, particle, model, params)
    total_energy = getEnergy(g) - thresh
    if rand() > 0.5 || :SP in model.meas || particle.type isa Positron
        other_part_energy = ion_type.loss(getEnergy(g)) #/params.energy_scale)*params.energy_scale*q
    else
        other_part_energy = total_energy - ion_type.loss(getEnergy(g)) #/params.energy_scale)*params.energy_scale*q
    end
    part_energy = total_energy - other_part_energy

    if part_energy < 0
        @info "part" part_energy #/params.energy_scale/q
        @info "other" other_part_energy #/params.energy_scale/q
    end

    return part_energy, other_part_energy
end

# function apply_attachment_weighting!(::Positron, particle)
#     particle.alive = false
#     particle
# end

function apply_attachment_weighting!(::Lepton, particle)
    particle.weight /= 2
    return particle
end

function attachmentColl!(particle, method)
    if method == :WEIGHTING
        apply_attachment_weighting!(particle.type, particle)
    elseif method == :EXPLICIT
        particle.alive = false
        particle.weight = 0.0
    else
        error("Unknown non-conservative method in ion_type")
    end
    return particle
end

function scatteringAngle(cs_no, cs, mass, g, particle, model, params)
    coschi = params.DCS[cs_no](getEnergy(g), rand())
    #' coschi = params.splat_gas.dcs[cs_no][particle.cbin](rand())::Float64
    #coschi = -1
    sinchi = sqrt(1 - coschi^2)
    phi = rand() * 2 * pi
    #' @debug """
    #' Scattering angle:
    #' $coschi $sinchi $phi"
    #' """
    return coschi, sinchi, phi
end

#@fastmath function shift(v_new::SVector{3, Quantity},v::SVector{3, Quantity})
@fastmath function shift(v_new, v)
    relvel_mag = norm(v)
    if relvel_mag == 0
        orig_costheta = 1.0
        orig_sintheta = 0.0
        orig_cosphi = 1.0
        orig_sinphi = 0.0
    else
        orig_costheta = v[3] / relvel_mag
        orig_sintheta = sqrt(1 - orig_costheta^2)
        if orig_sintheta == 0
            orig_cosphi = 1.0
            orig_sinphi = 0.0
        else
            orig_cosphi = v[1] / relvel_mag / orig_sintheta
            orig_sinphi = v[2] / relvel_mag / orig_sintheta
        end
    end

    Ry = @SMatrix [orig_costheta 0.0 orig_sintheta; 0.0 1.0 0.0; -orig_sintheta 0.0 orig_costheta]
    Rz = @SMatrix [orig_cosphi -orig_sinphi 0.0; orig_sinphi orig_cosphi 0.0; 0.0 0.0 1.0]

    return Rz * Ry * v_new
end

#### * Temperature
function findu(cs_ind, cs, mass, particle, model, params)
    #TODO: make generic function based on inputs
    w = sqrt(2 * k * model.T.val / mass)
    v = particle.v #*params.time_scale/params.length_scale
    he = uconvert(eV, electron_v2eV * (norm(v) + 5.0 * w)^2)
    le = uconvert(eV, max(1.0e-25u"kg*m^2/s^2", electron_v2eV * (norm(v) - 5.0 * w) .^ 2))

    #TODO: get actual CF from CS

    hcf = νₜ(params.gas, sqrt(he))
    lcf = νₜ(params.gas, sqrt(le))

    mcf = max(hcf, lcf)
    while true
        # usually it only goes through one while loop iteration
        #' u = (@SVector randn(3)) * (w / sqrt(2))
        u = SVector(randn(), randn(), randn()) * (w / sqrt(2))

        g = v .- u
        # @show v
        # @show u
        # @show g
        #' eg = uconvert(eV, electron_v2eV * dot(g, g))
        # @show eg
        ccf = νₜ(params.gas, sqrt(uconvert(eV, electron_v2eV * dot(g, g))))
        # @show lcf
        # @show hcf
        # @show ccf
        # @show mcf

        if ccf / mcf == 0.0
            error("Zero relative energy")
        end
        if rand() < ccf / mcf
            return g
        end
    end
    #' @debug "Neutral velocity: x = $(g[1]), y = $(g[2]), z = $(g[3])"
    return
end
