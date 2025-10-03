using Distributed
addprocs()
#' @everywhere using Pkg
#' @everywhere Pkg.activate("LightLeptonMC")

@info "Using $(nprocs()) cores."
using Test

using LightLeptonMC, LightLeptonTools
using DataFrames, DataFramesMeta, Statistics
function compareSwarm(A, B)
    for row in 1:nrow(A)
        for (i, t) in enumerate(A.types[row])
            @info """
            type:    $(t)
            Ref:     $(A.ys[row][i][1]|>upreferred)
            Current: $(B.ys[row][i][1]|>upreferred)
            """
            @test A.ys[row][i] ≈ B.ys[row][i] rtol = 1e-2
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
        sol = LightLeptonMC.mc(MaxwellModel(), swarmData, N=1_000_000, t=3e-5s)
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
        sol = LightLeptonMC.mc(MaxwellModel(), swarmData, N=1_000_000, t=3e-5s)
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
        sol = LightLeptonMC.mc(HardSphereModel(), swarmData, N=100_000, t=1e-5s)
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
        sol = LightLeptonMC.mc(HardSphereModel(), swarmData, N=100_000, t=1e-5s)
        compareSwarm(swarmData, sol)
    end
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
#'         sol = pt(MaxwellModel(), E=1.0, T=0.0, lmax=9, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "ReidRamp" begin
#'         @info "Inelastic (ReidRamp) Test"
#'         meas, swarm, params = mc(ReidRamp(), E=1.0, T=0.0, N=N, t=2e-7, meas=[:f, :t])
#'         sol = pt(ReidRamp(), E=1.0, T=0.0, lmax=9, N=5000)
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
#'         sol = pt(MaxwellModel(), E=1.0, T=T, lmax=9, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "Lucas-Saelee" begin
#'         @info "Ionisation (Lucas-Saelee) Test"
#'         T = 0.0
#'         meas, swarm, params = mc(LucasSaelee(), E=10.0, T=T, N=N, t=5e-7, meas=[:f, :t])
#'         sol = pt(LucasSaelee(), E=10.0, T=T, lmax=9, N=5000)
#'         compare(sol, swarm)
#'     end
#'     @testset "Ness-Robson" begin
#'         @info "Attachment (Ness-Robson) Test"
#'         T = 293.0
#'         meas, swarm, params = mc(PowerLawAttachment(), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(PowerLawAttachment(), E=0.4, T=T, lmax=9, N=5000)
#'         compare(sol, swarm)
#'         meas, swarm, params = mc(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, lmax=9, N=5000)
#'         compare(sol, swarm)
#'         meas, swarm, params = mc(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
#'         sol = pt(PowerLawAttachment(a=1 // 1000), E=0.4, T=T, lmax=9, N=5000)
#'         compare(sol, swarm)
#'     end
#' end
