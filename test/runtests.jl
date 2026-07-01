using Distributed
addprocs(20)
@everywhere using Pkg
@everywhere Pkg.activate("..")

@everywhere using Test

@everywhere using MCSim, Gases

N = 100000
tol = 0.01
function withinTol(y, ŷ, tol=0.05)
    @info "Values: " y ŷ
    @test y - y * tol <= ustrip(ŷ) <= y + y * tol
end

function compare(sol, swarm)
    withinTol(sol.ε̄, swarm.ε)
    withinTol(sol.W, swarm.bW)
    withinTol(sol.WF, swarm.fW)
    withinTol(sol.DT, swarm.bDt)
    withinTol(sol.DL, swarm.bDl)
    withinTol(sol.DTF, swarm.fDt)
    withinTol(sol.DLF, swarm.fDl)
    withinTol(sol.kio - sol.kat, swarm.R)
end

@testset "Swarm Benchmarks" begin
    @testset "Maxwell" begin
        @info "Elastic (Maxwell) Test"
        meas, swarm, params = mc(MaxwellModel(), E=1.0, T=0.0, N=N, t=3e-5, meas=[:f, :t])
        sol = pt(MaxwellModel(), E=1.0, T=0.0, lmax=9, N=5000)
        compare(sol, swarm)
    end
    @testset "ReidRamp" begin
        @info "Inelastic (ReidRamp) Test"
        meas, swarm, params = mc(ReidRamp(), E=1.0, T=0.0, N=N, t=2e-7, meas=[:f, :t])
        sol = pt(ReidRamp(), E=1.0, T=0.0, lmax=9, N=5000)
        compare(sol, swarm)
    end
    @testset "Anisotropic (ReidRamp)" begin
        @info "Anisotropic (ReidRamp) Test"
        meas, swarm, params = mc(ReidRampA(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
        withinTol(1.234, swarm.ε)
        withinTol(1.958e24, swarm.bDt)
        meas, swarm, params = mc(ReidRampB(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
        withinTol(1.229, swarm.ε)
        withinTol(1.898e24, swarm.bDt)
        meas, swarm, params = mc(ReidRampC(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
        withinTol(1.216, swarm.ε)
        withinTol(1.919e24, swarm.bDt)
        meas, swarm, params = mc(ReidRampD(), E=25.0, T=0.0, N=N, t=1e-7, meas=[:f, :t])
        withinTol(1.210, swarm.ε)
        withinTol(1.795e24, swarm.bDt)
    end
    @testset "Maxwell - Temperature" begin
        @info "Temperature (ReidRamp) Test"
        T = 298.0
        meas, swarm, params = mc(MaxwellModel(), E=1.0, T=T, N=N, meas=[:f, :t])
        sol = pt(MaxwellModel(), E=1.0, T=T, lmax=9, N=5000)
        compare(sol, swarm)
    end
    @testset "Lucas-Saelee" begin
        @info "Ionisation (Lucas-Saelee) Test"
        T = 0.0
        meas, swarm, params = mc(LucasSaelee(), E=10.0, T=T, N=N, t=5e-7, meas=[:f, :t])
        sol = pt(LucasSaelee(), E=10.0, T=T, lmax=9, N=5000)
        compare(sol, swarm)
    end
    @testset "Ness-Robson" begin
        @info "Attachment (Ness-Robson) Test"
        T = 293.0
        meas, swarm, params = mc(PowerLawAttachment(), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
        sol = pt(PowerLawAttachment(), E=0.4, T=T, lmax=9, N=5000)
        compare(sol, swarm)
        meas, swarm, params = mc(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
        sol = pt(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, lmax=9, N=5000)
        compare(sol, swarm)
        meas, swarm, params = mc(PowerLawAttachment(a=1 // 10000), E=0.4, T=T, N=N, t=3e-5, meas=[:f, :t])
        sol = pt(PowerLawAttachment(a=1 // 1000), E=0.4, T=T, lmax=9, N=5000)
        compare(sol, swarm)
    end
end
