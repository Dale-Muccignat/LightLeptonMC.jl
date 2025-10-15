function initParams(input_gas=MaxwellModel(), model::MCMODEL=MCMODEL())

    #' cbins = [0.0, (10.0 .^ (-8:0.1:4))...] .* eV
    cbins = [0.0, (10.0 .^ (-8:1:-2))..., (10.0 .^ (-1.9:0.1:4))...] .* eV

    tbins, ebins, rbins, χbins, zbins, ybins, xbins = generateBins(model)
    n0 = buildDensity(model.n0, rbins, zbins)
    a = buildField(model.En0)

    total_freq = upreferred.(buildFreq(input_gas, cbins, model))
    # dcs = buildDCS(input_gas, cbins, model)

    params = MCPARAMS(
        ebins=uconvert.(J̃, ebins),
        tbins=upreferred.(tbins),
        zbins=upreferred.(zbins),
        ybins=upreferred.(ybins),
        xbins=upreferred.(xbins),
        cbins=uconvert.(J̃, cbins),
        rbins=upreferred.(rbins),
        χbins=upreferred.(χbins),
        gas=input_gas,
        n0=upreferred.(n0),
        a=upreferred.(a),
        total_freq=total_freq,
    )

    return params
end

function initMeas(model, params::MCPARAMS=MCPARAMS())
    MEAS(model.type, model, params)
end

#### DCS ####

function buildDCS(σ::DifferentialCrossSection, cbins)
    # Create a 2d matrix interpolation that takes ε and R as an input and returns χ
    Rs = range(0, 1, 100)
    χs = range(0, π, 100)
    out = Matrix{typeof{1.0u"rad"}}(undef, length(cbins), length(Rs))

    for (k, (g, _, _)) in enumerate(gas)
        # create a minimising function
        for (i, ε) in enumerate(cbins)
            total = quadgk((χ) -> g(ε, χ) * sin(χ))[1]
            for (j, r) in enumerate(Rs)
                L(χ) = abs(r - quadgk((χ) -> g(ε, χ) * sin(χ) / total, 0, χ)[1])
                out[i, j, k] = cos(Optim.minimizer(optimize(L, 0, 1)))
            end
        end
    end
    linear_interpolation((cbins, Rs), out)
end

function buildDCS(σ::AbstractCrossSection, cbins)
    # Create a 2d matrix interpolation that takes ε and R as an input and returns χ
    (R) -> 1 - 2 * R
end

#### Frequency ####

function buildFreq(gas::LLCS, cbins, model)
    total_freq = zeros(length(cbins)) .* sqrt(eV) .* Å^2
    en_spacing = cbins[end] - cbins[end-1]
    inv_total_weight(ε) = -νₜ(gas, ε) - sum(CollisionFrequency(σ)(sqrt(ε)) for σ in gas.attachment)
    inv_total(ε) = -νₜ(gas, ε)

    if model.ion_type.method === :WEIGHTING
        for i in 1:length(total_freq)-1
            total_freq[i] = -Optim.minimum(optimize(inv_total_weight, ustrip(sqrt(eV), sqrt(cbins[i])), ustrip(sqrt(eV), sqrt(cbins[i+1])))) * sqrt(eV) * Å^2
        end
        total_freq[end] = -Optim.minimum(optimize(inv_total_weight, ustrip(sqrt(eV), sqrt(cbins[end])), ustrip(sqrt(eV), sqrt(cbins[end] + en_spacing)))) * sqrt(eV) * Å^2
    else
        for i in 1:length(total_freq)-1
            total_freq[i] = -Optim.minimum(optimize(inv_total, ustrip(sqrt(eV), sqrt(cbins[i])), ustrip(sqrt(eV), sqrt(cbins[i+1])))) * sqrt(eV) * Å^2
        end
        total_freq[end] = -Optim.minimum(optimize(inv_total, ustrip(sqrt(eV), sqrt(cbins[end])), ustrip(sqrt(eV), sqrt(cbins[end] + en_spacing)))) * sqrt(eV) * Å^2
    end

    return total_freq
end

#### * Bins
# This tries a best guess at what parameters will be needed but is in no way good. Once params are initialized you will need to customize the ones you want as you need

function generateBins(model)

    m = string.(model.meas)

    tbins = 1 ∈ occursin.("t", m) ? buildTBins(model.type, model) : nothing
    ebins = 1 ∈ occursin.("ε", m) ? buildEnBins(model.type, model) : nothing
    zbins = 1 ∈ occursin.("z", m) ? buildZBins(model.type, model) : nothing
    ybins = 1 ∈ occursin.("y", m) ? buildYBins(model.type, model) : nothing
    xbins = 1 ∈ occursin.("x", m) ? buildXBins(model.type, model) : nothing
    rbins = 1 ∈ occursin.("r", m) ? buildRBins(model.type, model, model.n0) : nothing
    χbins = 1 ∈ occursin.("χ", m) ? buildχBins(model.type, model) : nothing

    return tbins, ebins, rbins, χbins, zbins, ybins, xbins

end

# TODO: Combine these into one function per model

function buildRBins(type::SWARM, model, n0::UniformDensity)
    return LinRange(0m, 1e-6m, 10)
end

function buildEnBins(type::SWARM, model)
    return 10 .^ LinRange(-3, 2, 1000) .* eV
end

function buildTBins(type::SWARM, model)
    return LinRange(0s, type.t_f, 200)
end

function buildχBins(type::SWARM, model)
    return LinRange(0, π, 90) .* u"rad"
end

function buildZBins(type::SWARM, model)
    return LinRange(-10mm, 10mm, 200)
end

function buildYBins(type::SWARM, model)
    return LinRange(-7e-7, 7e-7, 200) .* m
end

function buildXBins(type::SWARM, model)
    return LinRange(-7e-7, 7e-7, 200) .* m
end

function buildRBins(type::CELL, model, n0)
    return LinRange(0, type.r, type.bin_no) .* m^2
end

function buildEnBins(type::CELL, model)
    #TODO: Make specific to each velocity type
    #return collect(model.type.thresh:model.type.Δε_detector:(model.v_0.ε+2*model.v_0.Δε_source))
    # return 10 .^LinRange(-3,3,600) .*eV
    return collect(10.0 .^ (-3:0.1:6)) .* eV
end

function buildTBins(type::CELL, model)
    # return LinRange(0s,type.t_f,10000)
    return [-Inf * s, (0s:1e-9s:type.t_f)...]
end

function buildχBins(type::CELL, model)
    return LinRange(0, π, 100) .* u"rad"
end

function buildZBins(type::CELL, model)
    # return [-10m, LinRange(0. *m,type.z,type.z_n)...]
    return LinRange(0.0 * m, type.z, type.z_n) |> collect
end

function buildYBins(type::CELL, model)
    return LinRange(-type.y, type.y, 200)
end

function buildXBins(type::CELL, model)
    return LinRange(-type.x, type.x, 200)
end

function buildRBins(type::LμJ, model, n0::LμJDensity)
    return [0.0 * m^2, (type.jet_width / 2)^2, 0.02^2 * m^2]
end

function buildRBins(type::LμJ, model, n0::LμJDensityGradient)
    return adaptiveMeshlog(n0.func, 0.0 * m, 0.02 * m, n0.δtol) .^ 2
end

function buildEnBins(type::LμJ, model)
    return collect(model.type.thresh:0.3eV:(model.v_0.ε+4*model.v_0.Δε_source))
end

function buildTBins(type::LμJ, model)
    return LinRange(3e-13, 320e-13, 10) .* s
end

function buildχBins(type::LμJ, model)
    return LinRange(0, π, 100) .* u"rad"
end

function buildZBins(type::LμJ, model)
    return LinRange(-7.5e-6, -6.5e-6, 200) .* m
end

function buildYBins(type::LμJ, model)
    return LinRange(-7e-7, 7e-7, 200) .* m
end

function buildXBins(type::LμJ, model)
    return LinRange(-7e-7, 7e-7, 200) .* m
end

#### Density ####

function buildDensity(n0::UniformDensity, rbins, zbins)
    n0.val
end

function buildDensity(n0::LμJDensityGradient, rbins, zbins)
    n0.func.(sqrt.(rbins))
end

function buildDensity(n0::LμJDensity, rbins, zbins)
    [n0.val, n0.background, 0.0 * m^-3]
end

function buildDensity(n0::CellDensity, rbins, zbins)
    n0.val .* ones(length(zbins))
end

#### Electric Field ####

function buildField(En0::UniformField)
    SVector(0.0 * m^4 * s^-2, 0.0 * m^4 * s^-2, upreferred((q * En0.val / me)))
end

function buildField(En0::ZeroField)
    SVector(0.0 * m^4 / s^-2, 0.0 * m^4 / s^-2, 0.0 * m^4 / s^-2)
end


#### * PARTICLE

function PART(model, params)

    particle = PART()
    particle.type = model.particle_type
    particle.r = upreferred.(initPartr(model.r_0, params))#.*params.length_scale
    particle.v = upreferred.(initPartv(model.v_0, params))#.*params.length_scale ./params.time_scale
    particle.t = upreferred.(initPartt(model.t_0, params))
    particle.cbin = findMinInd(params.cbins, getEnergy(particle.v))
    if params.ebins != nothing
        particle.ebin = findMinInd(params.ebins, getEnergy(particle.v))
    end
    if params.tbins != nothing
        particle.tbin = findMinInd(params.tbins, particle.t)
    end

    if params.rbins != nothing
        particle.rbin = findMinInd(params.rbins, dot(particle.r, particle.r))
    end

    if params.zbins != nothing
        particle.zbin = findMinInd(params.zbins, particle.r[3])
    end

    if params.xbins != nothing
        particle.xbin = findMinInd(params.xbins, particle.r[1])
    end

    if params.ybins != nothing
        particle.ybin = findMinInd(params.ybins, particle.r[2])
    end


    return particle

end



function initPartr(r_0::ZeroSpatial, params)
    r = SArray{Tuple{3},Float64,1,3}(0.0, 0.0, 0.0)
    return r * m
end
function initPartr(r_0::InitSpatial, params)
    r = r_0.axis .* r_0.val
end
function initPartr(r_0::FunctionalSpatial, params)
    r_0.getAxis() .* r_0.getPosition(rand())
end
function initPartr(r_0::GaussianLμJSpatial, params)
    # SD = r_0.Δy_source/(2*sqrt(2*log(2)))

    return r_0.Δy_source / (2 * sqrt(2 * log(2))) * randn() .* r_0.gaussian_axis .+ r_0.val .* r_0.axis
    # val = Inf*m .* r_0.gaussian_axis .+ r_0.val .* r_0.axis
    # while -7.5μm > val[2] || val[2] > 7.5μm 
    #     R=randn()
    #     # while val < -SD*R
    #         # R=randn()
    #     # end
    #     val= SD*R .* r_0.gaussian_axis .+ r_0.val .* r_0.axis
    # end
    # val
end

function initPartt(t_0::InitTemporal, params)
    t_0.val
end

function initPartt(t_0::GaussianTemporal, params)
    # electron_v2eV = 0.5 * me
    # coschi = 1-2*rand()
    # sinchi = sqrt(1-coschi^2)
    # phi = rand()*2*pi
    #
    return t_0.Δt_source / (2 * sqrt(2 * log(2))) * randn() .+ t_0.val
end

function initPartv(v_0::ZeroVelocity, params)
    r = SArray{Tuple{3},Float64,1,3}(0.0, 0.0, 0.0)
    return r * m / s
end

function initPartv(v_0::FunctionalEnergy, params)
    #TODO: Make all these static arrays
    return v_0.getAxis() .* getVelocity(v_0.getVelocity(rand()))
end

function initPartv(v_0::InitVelocity, params)
    #TODO: Make all these static arrays
    v = v_0.axis .* getVelocity(v_0.ε)
    return v
end

function initPartv(v_0::GaussianDirectionVelocity, params)
    electron_v2eV = 0.5 * me
    coschi = 1 - 2 * rand()
    sinchi = sqrt(1 - coschi^2)
    phi = rand() * 2 * pi

    sqrt(v_0.ε / electron_v2eV) * SArray{Tuple{3},Float64,1,3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
end


function initPartv(v_0::GaussianEnergyVelocity, params)
    SD = v_0.Δε_source / (2 * sqrt(2 * log(2)))
    electron_v2eV = 0.5 * me

    ea = v_0.ε

    R = randn()
    while ea < -SD * R
        R = randn()
    end
    new_ea = ea + SD * R

    sqrt.(new_ea / electron_v2eV) * v_0.axis
end

