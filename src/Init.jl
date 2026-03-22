function initParams(input_gas = MaxwellModel(), model::MCMODEL = MCMODEL())

    #' cbins = [0.0, (10.0 .^ (-8:0.1:4))...] .* eV
    cbins = [0.0, (10.0 .^ (-8:1:-2))..., (10.0 .^ (-1.9:0.1:4))...] .* eV

    tbins, ebins, rbins, χbins, zbins, ybins, xbins = generateBins(model)
    n0 = buildDensity(model.n0, rbins, zbins)
    a = buildField(model.En0)

    total_freq = upreferred.(buildFreq(input_gas, cbins, model))
    dcs = [buildDCS(g, cbins) for (g, _, _) in input_gas]

    params = MCPARAMS(
        ebins = uconvert.(J̃, ebins),
        tbins = upreferred.(tbins),
        zbins = upreferred.(zbins),
        ybins = upreferred.(ybins),
        xbins = upreferred.(xbins),
        cbins = uconvert.(J̃, cbins),
        rbins = upreferred.(rbins),
        χbins = upreferred.(χbins),
        gas = input_gas,
        n0 = upreferred.(n0),
        a = upreferred.(a),
        total_freq = total_freq,
        DCS = dcs,
    )

    return params
end

function initMeas(model, params::MCPARAMS = MCPARAMS())
    return MEAS(model.type, model, params)
end

#### DCS ####

@inbounds function buildDCS(σ::DifferentialCrossSection, cbins = logrange(1.0e-2, 1.0e2, 41) .* eV)
    # Create a 2d matrix interpolation that takes ε and R as an input and returns χ
    Rs = range(0, 1, 25)
    χs = range(0.0 * u"rad", π * u"rad", 25)
    out = Matrix{Float64}(undef, length(cbins), length(Rs))
    c1, c2 = χs[1], χs[end]
    spacing = χs[2] - χs[1]
    # create a minimising function
    last = 0.0
    for i in eachindex(cbins)
        function intχ(χ)
            σ(cbins[i], χ) * sin(χ)
        end
        total = quadgk(intχ, c1, c2)[1]::typeof(1.0u"Å^2*rad")
        for j in eachindex(Rs)
            function L(χ)
                partial = quadgk(intχ, c1, χ * u"rad")[1]::typeof(1.0u"Å^2*rad")
                abs(Rs[j] - partial / total)
            end
            result = optimize(L, last, π)
            last = Optim.minimizer(result)
            out[i, j] = cos(last)
        end
        last = 0.0
    end
    linear_interpolation((cbins, Rs), out, extrapolation_bc = Interpolations.Flat())
end

function buildDCS(σ::T, cbins = logrange(1.0e-2, 1.0e2, 41) .* eV) where {T <: AbstractCrossSection}
    # Create a 2d matrix interpolation that takes ε and R as an input and returns χ
    return (ε, R) -> 1 - 2 * R
end

#### Liquid ####

@inbounds function buildLiquid(σ::AbstractCrossSection, S, cbins)
    # Takes an elastic cross section, σ, and an angle dependent static structure factor S
    # Returns three cross sections, both, energy, momentum
    function int(χ, ε)
        0.5 * S(2 * sqrt(2 * me * ε) / ħ * sqrt((1 - cos(χ)) / 2)) * (1 - cos(χ)) * sin(χ)
    end
    σ_b = Vector{uÅ2}(undef, length(cbins))
    σ_e = Vector{uÅ2}(undef, length(cbins))
    σ_m = Vector{uÅ2}(undef, length(cbins))
    for i in eachindex(cbins)
        e = cbins[i]
        res = quadgk((x) -> int(x, e), 0, pi)[1]::Float64
        @show res
        if res >= 1
            σ_b[i] = σ(e)
            σ_e[i] = 0.0u"Å^2"
            σ_m[i] = σ(e) * (res - 1)
        else
            σ_b[i] = σ(e) * res
            σ_e[i] = σ(e) * (1 - res)
            σ_m[i] = 0.0u"Å^2"
        end
    end
    return CrossSection(collect(cbins), σ_b; type = :ELASTIC), CrossSection(collect(cbins), σ_e; type = :ELASTIC_ENERGY), CrossSection(collect(cbins), σ_m; type = :ELASTIC_MOMENTUM)
end
#### Frequency ####

function buildFreq(gas::LLCS, cbins, model)
    total_freq = zeros(length(cbins)) .* sqrt(eV) .* Å^2
    en_spacing = cbins[end] - cbins[end - 1]
    inv_total_weight(ε) = -νₜ(gas, ε) - sum(CollisionFrequency(σ)(sqrt(ε)) for σ in gas.attachment; init = 0.0)
    inv_total(ε) = -νₜ(gas, ε)

    if model.ion_type.method === :WEIGHTING
        for i in 1:(length(total_freq) - 1)
            total_freq[i] = -Optim.minimum(optimize(inv_total_weight, ustrip(sqrt(eV), sqrt(cbins[i])), ustrip(sqrt(eV), sqrt(cbins[i + 1])))) * sqrt(eV) * Å^2
        end
        total_freq[end] = -Optim.minimum(optimize(inv_total_weight, ustrip(sqrt(eV), sqrt(cbins[end])), ustrip(sqrt(eV), sqrt(cbins[end] + en_spacing)))) * sqrt(eV) * Å^2
    else
        for i in 1:(length(total_freq) - 1)
            total_freq[i] = -Optim.minimum(optimize(inv_total, ustrip(sqrt(eV), sqrt(cbins[i])), ustrip(sqrt(eV), sqrt(cbins[i + 1])))) * sqrt(eV) * Å^2
        end
        total_freq[end] = -Optim.minimum(optimize(inv_total, ustrip(sqrt(eV), sqrt(cbins[end])), ustrip(sqrt(eV), sqrt(cbins[end] + en_spacing)))) * sqrt(eV) * Å^2
    end

    return total_freq
end

#### * Bins
# This tries a best guess at what parameters will be needed but is in no way good. Once params are initialized you will need to customize the ones you want as you need

function generateBins(model)

    m_bins = string.(model.meas)

    tbins = 1 ∈ occursin.("t", m_bins) ? buildTBins(model.type, model) : LinRange(0.0s, 0.0s, 0)
    ebins = 1 ∈ occursin.("ε", m_bins) ? buildEnBins(model.type, model) : LinRange(0.0eV, 0.0eV, 0)
    zbins = 1 ∈ occursin.("z", m_bins) ? buildZBins(model.type, model) : LinRange(0.0m, 0.0m, 0)
    ybins = 1 ∈ occursin.("y", m_bins) ? buildYBins(model.type, model) : LinRange(0.0m, 0.0m, 0)
    xbins = 1 ∈ occursin.("x", m_bins) ? buildXBins(model.type, model) : LinRange(0.0m, 0.0m, 0)
    rbins = 1 ∈ occursin.("r", m_bins) ? buildRBins(model.type, model, model.n0) : LinRange(0.0m^2, 0.0m^2, 0)
    χbins = 1 ∈ occursin.("χ", m_bins) ? buildχBins(model.type, model) : LinRange(0.0 * u"rad", 0.0 * u"rad", 0)

    return tbins, ebins, rbins, χbins, zbins, ybins, xbins

end

# TODO: Combine these into one function per model

function buildRBins(type::SWARM, model, n0::UniformDensity)
    return LinRange(0m^2, 1.0e-6m^2, 10)
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
    return LinRange(-7.0e-7, 7.0e-7, 200) .* m
end

function buildXBins(type::SWARM, model)
    return LinRange(-7.0e-7, 7.0e-7, 200) .* m
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
    return [-Inf * s, (0s:1.0e-9s:type.t_f)...]
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
    return collect(model.type.thresh:model.type.Δε_detector:(model.v_0.ε + 4 * model.v_0.Δε_source))
end

function buildTBins(type::LμJ, model)
    return LinRange(3.0e-13, 320.0e-13, 10) .* s
end

function buildχBins(type::LμJ, model)
    return LinRange(0, π, 100) .* u"rad"
end

function buildZBins(type::LμJ, model)
    return LinRange(-7.5e-6, -6.5e-6, 200) .* m
end

function buildYBins(type::LμJ, model)
    return LinRange(-7.0e-7, 7.0e-7, 200) .* m
end

function buildXBins(type::LμJ, model)
    return LinRange(-7.0e-7, 7.0e-7, 200) .* m
end

#### Density ####

function buildDensity(n0::UniformDensity, rbins, zbins)
    return n0.val
end

function buildDensity(n0::LμJDensityGradient, rbins, zbins)
    return n0.func.(sqrt.(rbins))
end

function buildDensity(n0::LμJDensity, rbins, zbins)
    return [n0.val, n0.background, 0.0 * m^-3]
end

function buildDensity(n0::CellDensity, rbins, zbins)
    return n0.val .* ones(length(zbins))
end

#### Electric Field ####

function buildField(En0::UniformField)
    return SVector(0.0 * m^4 * s^-2, 0.0 * m^4 * s^-2, upreferred((q * En0.val / me)))
end

function buildField(En0::ZeroField)
    return SVector(0.0 * m^4 / s^-2, 0.0 * m^4 / s^-2, 0.0 * m^4 / s^-2)
end


#### * PARTICLE

function PART(model, params)

    particle = PART()
    particle.type = model.particle_type
    particle.r = upreferred.(initPartr(model.r_0, params)) #.*params.length_scale
    particle.v = upreferred.(initPartv(model.v_0, params)) #.*params.length_scale ./params.time_scale
    particle.t = upreferred.(initPartt(model.t_0, params))
    particle.cbin = findMinInd(params.cbins, getEnergy(particle.v))
    if params.ebins |> !isempty
        particle.ebin = findMinInd(params.ebins, getEnergy(particle.v))
    end
    if params.tbins |> !isempty
        particle.tbin = findMinInd(params.tbins, particle.t)
    end

    if params.rbins |> !isempty
        particle.rbin = findMinInd(params.rbins, dot(particle.r, particle.r))
    end

    if params.zbins |> !isempty
        particle.zbin = findMinInd(params.zbins, particle.r[3])
    end

    if params.xbins |> !isempty
        particle.xbin = findMinInd(params.xbins, particle.r[1])
    end

    if params.ybins |> !isempty
        particle.ybin = findMinInd(params.ybins, particle.r[2])
    end


    return particle

end


function initPartr(r_0::ZeroSpatial, params)
    r = SArray{Tuple{3}, Float64, 1, 3}(0.0, 0.0, 0.0)
    return r * m
end
function initPartr(r_0::InitSpatial, params)
    return r = r_0.axis .* r_0.val
end
function initPartr(r_0::FunctionalSpatial, params)
    return r_0.getAxis() .* r_0.getPosition(rand())
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
    return t_0.val
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
    r = SArray{Tuple{3}, Float64, 1, 3}(0.0, 0.0, 0.0)
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

    return sqrt(v_0.ε / electron_v2eV) * SArray{Tuple{3}, Float64, 1, 3}(cos(phi) * sinchi, sin(phi) * sinchi, coschi)
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

    return sqrt.(new_ea / electron_v2eV) * v_0.axis
end
