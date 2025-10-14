using Distributed
addprocs()
#' @everywhere using Pkg
#' @everywhere Pkg.activate("LightLeptonMC")

@info "Using $(nprocs()) cores."
using Test

using LightLeptonMC, LightLeptonTools
using DataFrames, DataFramesMeta, Statistics, Measurements

function compareSwarm(A, B)
    for row in 1:nrow(A)
        for (i, t) in enumerate(A.types[row])
            @info """
            type:    $(t)
            Ref:     $(A.ys[row][i][1])
            Current: $(B.ys[row][i][1])
            """
            @test A.ys[row][i][1] ≈ Measurements.value(B.ys[row][i][1])  rtol = 1e-2
        end
    end
end

@testset "Maxwell Model" begin
    # Write your tests here

    @testset "Maxwell Model - 0K" begin
        # Reference: Boyle, 2016, PhD Thesis.
        swarmData = DataFrame(
            x=[[1.0Td]],
            T=0.0K,
            P=0.0kPa,
            types=[[:ε, :W, :nDT, :nDL]],
            ys=[[
                [0.72916eV],
                [5931.0m / s],
                [2.8831e24 / s / m],
                [2.8831e24 / s / m],
            ]]
        )
        sol = LightLeptonMC.mc(MaxwellModel(), swarmData, N=500_000, t=3e-5s, S=0.1eV)
        compareSwarm(swarmData, sol)
    end

    @testset "Maxwell Model - 293K" begin
        # Reference: Boyle, 2016, PhD Thesis.
        swarmData = DataFrame(
            x=[[1.0Td]],
            T=293.0K,
            P=0.0kPa,
            types=[[:ε, :W, :nDT, :nDL]],
            ys=[[
                [0.76703eV],
                [5931.0m / s],
                [3.0328e24 / s / m],
                [3.0328e24 / s / m],
            ]]
        )
        sol = LightLeptonMC.mc(MaxwellModel(), swarmData, N=500_000, t=3e-5s, S=0.1eV)
        compareSwarm(swarmData, sol)
    end
end

@testset "HardSphere Model" begin
    # Reference: Boyle, 2016, PhD Thesis.
    @testset "HardSphere Model - 0K" begin
        swarmData = DataFrame(
            x=[[1.0Td]],
            T=0.0K,
            P=0.0kPa,
            types=[[:ε, :W, :nDT, :nDL]],
            ys=[[
                [0.72932eV],
                [5759.5m / s],
                [3.202e24 / s / m],
                [1.573e24 / s / m],
            ]]
        )
        sol = LightLeptonMC.mc(HardSphereModel(), swarmData, N=500_000, t=2e-5s)
        compareSwarm(swarmData, sol)
    end

    @testset "HardSphere Model - 293K" begin
        # Reference: Boyle, 2016, PhD Thesis.
        swarmData = DataFrame(
            x=[[1.0Td]],
            T=293.0K,
            P=0.0kPa,
            types=[[:ε, :W, :nDT, :nDL]],
            ys=[[
                [0.75057eV],
                [5639.6m / s],
                [3.247e24 / s / m],
                [1.640e24 / s / m],
            ]]
        )
        sol = LightLeptonMC.mc(HardSphereModel(), swarmData, N=500_000, t=2e-5s)
        compareSwarm(swarmData, sol)
    end
end

@testset "Power Law Attachment - 293K" begin
    # Reference: Boyle, 2016, PhD Thesis.
    swarmData = DataFrame(
        x=[[0.4Td]],
        T=293.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL, :WF, :nDTF, :nDLF, :Keff]],
        ys=[[
            [0.1364eV],
            [937.5m / s],
            [2.414e23 / s / m],
            [1.73e23 / s / m],
            [1186m / s],
            [2.697e23 / s / m],
            [2.697e23 / s / m],
            [-8.092e-19 * m^3 / s],
        ]]
    )
    g = PowerLawAttachment(a=0.001, p=0.5)
    sol = LightLeptonMC.mc(g, swarmData, N=1_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation(), S=0.1eV)
    sol = LightLeptonMC.mc(g, swarmData, N=500_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation(), S=0.1eV)
    compareSwarm(swarmData, sol)

    # # Reference: Boyle, 2016, PhD Thesis.
    # NOTE: This is too intensive
    #
    # swarmData = DataFrame(
    #     x=[[0.4Td]],
    #     T=293.0K,
    #     P=0.0kPa,
    #     types=[[:ε, :W, :nDT, :nDL, :WF, :nDTF, :nDLF, :Keff]],
    #     ys=[[
    #         [0.0847eV],
    #         [448m / s],
    #         [1.154e23 / s / m],
    #         [0.4072e23 / s / m],
    #         [1186m / s],
    #         [1.675e23 / s / m],
    #         [1.675e23 / s / m],
    #         [-5.026e-18 * m^3 / s],
    #     ]]
    # )
    # sol = LightLeptonMC.mc(PowerLawAttachment(a=0.01, p=0.5), swarmData, N=100_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation())
    # compareSwarm(swarmData, sol)

    #' @testset "Power Law - 293K, p=-1/2, a=0.01" begin
    # Reference: Boyle, 2016, PhD Thesis.
    swarmData = DataFrame(
        x=[[0.4Td]],
        T=293.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL, :WF, :nDTF, :nDLF, :Keff]],
        ys=[[
            [0.1545 * eV],
            [1186 * m / s],
            [3.055e23 / s / m],
            [3.055e23 / s / m],
            [1186 * m / s],
            [3.055e23 / s / m],
            [3.055e23 / s / m],
            [-5.931e-20 * m^3 / s],
        ]]
    )
    g = PowerLawAttachment(a=0.00001, p=-0.5)
    sol = LightLeptonMC.mc(g, swarmData, N=1_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation(), S=0.1eV)
    sol = LightLeptonMC.mc(g, swarmData, N=500_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation(), S=0.1eV)
    compareSwarm(swarmData, sol)

    #' @testset "Power Law - 293K, p=-1, a=0.00001" begin
    # Reference: Boyle, 2016, PhD Thesis.
    # NOTE: The effective attachment rate has a type-o in the thesis, powers of 10 should be positive
    swarmData = DataFrame(
        x=[[0.4Td]],
        T=293.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL, :WF, :nDTF, :nDLF, :Keff]],
        ys=[[
            [0.1572 * eV],
            [1227 * m / s],
            [3.159e23 / s / m],
            [3.199e23 / s / m],
            [1186 * m / s],
            [3.107e23 / s / m],
            [3.107e23 / s / m],
            [-2.054e-19 * m^3 / s],
        ]]
    )
    g = PowerLawAttachment(a=0.00001, p=-1)
    sol = LightLeptonMC.mc(g, swarmData, N=1_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation(), S=0.1eV)
    sol = LightLeptonMC.mc(g, swarmData, N=1_000_000, t=2e-5s, meas=[:f, :t], ion_type=LightLeptonMC.UniformIonisation(), S=0.1eV)
    compareSwarm(swarmData, sol)
end

@testset "Reid's Ramp" begin
    # Reference: Boyle, 2016, PhD Thesis.
    swarmData = DataFrame(
        x=[[1Td]],
        T=0.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL]],
        ys=[[
            [0.1015 * eV],
            [12720 * m / s],
            [9.75e23 / s / m],
            [7.591e23 / s / m],
        ]]
    )
    sol = LightLeptonMC.mc(ReidRamp(), swarmData, N=500_000, t=2e-7s)
    compareSwarm(swarmData, sol)

    # Reference: Boyle, 2016, PhD Thesis.
    swarmData = DataFrame(
        x=[[24Td]],
        T=0.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL]],
        ys=[[
            [0.4080 * eV],
            [88870 * m / s],
            [11.337e23 / s / m],
            [4.606e23 / s / m],
        ]]
    )
    sol = LightLeptonMC.mc(ReidRamp(), swarmData, N=500_000, t=2e-7s)
    compareSwarm(swarmData, sol)
end


@testset "Lucas and Saelee Ionisation" begin
    # Reference: Boyle, 2016, PhD Thesis.
    # NOTE: There is a type in the thesis. The α values are off by a factor of 100
    # F=1
    swarmData = DataFrame(
        x=[[10Td]],
        T=0.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL, :Keff]],
        ys=[[
            [4.969 * eV],
            [94740 * m / s],
            [27.23e24 / s / m],
            [29.33e24 / s / m],
            [2.42e-17 * m^3 / s],
        ]]
    )
    sol = LightLeptonMC.mc(LucasSaelee(F=1), swarmData, N=500_000, t=5e-7s)
    compareSwarm(swarmData, sol)

    # Reference: Boyle, 2016, PhD Thesis.
    # F=0.5
    swarmData = DataFrame(
        x=[[10Td]],
        T=0.0K,
        P=0.0kPa,
        types=[[:ε, :W, :nDT, :nDL, :Keff]],
        ys=[[
            [5.224 * eV],
            [85930 * m / s],
            [27.26e24 / s / m],
            [28.65e24 / s / m],
            [1.33e-17 * m^3 / s],
        ]]
    )
    sol = LightLeptonMC.mc(LucasSaelee(F=0.5), swarmData, N=500_000, t=5e-7s)
    compareSwarm(swarmData, sol)
end
#' N = 100000
#' tol = 0.01
#' function withinTol(y, ŷ, tol=0.05)
#'     @info "Values: " y ŷ
#'     @test y - y * tol <= ustrip(ŷ) <= y + y * tol
#' end

#' function compare(sol, swarm)
#'     withinTol(sol.ε̄, swarm.ε)
#'     withinTol(sol.W, swarm.bW)
#'     withinTol(sol.WF, swarm.fW)
#'     withinTol(sol.DT, swarm.bDt)
#'     withinTol(sol.DL, swarm.bDl)
#'     withinTol(sol.DTF, swarm.fDt)
#'     withinTol(sol.DLF, swarm.fDl)
#'     withinTol(sol.kio - sol.kat, swarm.R)
#' end

#' @testset "Swarm Benchmarks" begin
#'     @testset "Maxwell" begin
#'         @info "Elastic (Maxwell) Test"
#'         meas, swarm, params = mc(MaxwellModel(), E=1.0, T=0.0, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(MaxwellModel(), E=1.0, T=0.0, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "ReidRamp" begin
#'         @info "Inelastic (ReidRamp) Test"
#'         meas, swarm, params = mc(ReidRamp(), E=1.0, T=0.0, N=N, t=2e-7, meas=[:f, :t])
#'         sol = pt(ReidRamp(), E=1.0, T=0.0, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "Anisotropic (ReidRamp)" begin
#'         @info "Anisotropic (ReidRamp) Test"
#'         meas, swarm, params = mc(ReidRampA(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
#'         withinTol(1.234, swarm.ε)
#'         withinTol(1.958e24, swarm.bDt)
#'         meas, swarm, params = mc(ReidRampB(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
#'         withinTol(1.229, swarm.ε)
#'         withinTol(1.898e24, swarm.bDt)
#'         meas, swarm, params = mc(ReidRampC(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
#'         withinTol(1.216, swarm.ε)
#'         withinTol(1.919e24, swarm.bDt)
#'         meas, swarm, params = mc(ReidRampD(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
#'         withinTol(1.210, swarm.ε)
#'         withinTol(1.795e24, swarm.bDt)
#'     end
#'     @testset "Maxwell - Temperature" begin
#'         @info "Temperature (ReidRamp) Test"
#'         T = 298.0
#'         meas, swarm, params = mc(MaxwellModel(), E=1.0, T=T, N=N, meas=[:f, :t])
#'         sol = pt(MaxwellModel(), E=1.0, T=T, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "Lucas-Saelee" begin
#'         @info "Ionisation (Lucas-Saelee) Test"
#'         T = 0.0
#'         meas, swarm, params = mc(LucasSaelee(), E=10.0, T=T, N=N, t=5e-7, meas=[:f, :t])
#'         sol = pt(LucasSaelee(), E=10.0, T=T, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "Ness-Robson" begin
#'         @info "Attachment (Ness-Robson) Test"
#'         T = 293.0
#'         meas, swarm, params = mc(PowerLawAttachment(), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(PowerLawAttachment(), E=0.4, T=T, N=5000)
#'         compare(sol, swarm)
#'         meas, swarm, params = mc(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=5000)
#'         compare(sol, swarm)
#'         meas, swarm, params = mc(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(PowerLawAttachment(a=1 // 1000), E=0.4, T=T, N=5000)
#'         compare(sol, swarm)
#'     end
#' end
