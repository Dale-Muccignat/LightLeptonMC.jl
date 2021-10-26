using Distributed
addprocs(31)
using Test


using MCSim

N = 200000
tol = 0.01
function withinTol(y,ŷ,tol=0.01)
    @test y - y*tol <= ŷ<= y+y*tol
end
@testset "MCSim.jl" begin
    @testset "Maxwell" begin
        @info "Elastic (Maxwell) Test"
        meas, swarm, params  = mc(MaxwellModel(),E=1.,T=0.,N=N,t=3e-5)
        sol = pt(MaxwellModel(),E=1.,T=0.,lmax=9,N=5000)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "ReidRamp" begin
        @info "Inelastic (ReidRamp) Test"
        meas, swarm, params  = mc(ReidRamp(),E=1.,T=0.,N=N,t=2e-7)
        sol = pt(ReidRamp(),E=1.,T=0.,lmax=9,N=5000)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "Maxwell - Temperature" begin
        @info "Temperature (ReidRamp) Test"
        T = 298.
        meas, swarm, params  = mc(MaxwellModel(),E=1.,T=T,N=N)
        sol = pt(MaxwellModel(),E=1.,T=T,lmax=9,N=5000)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "Lucas-Saelee" begin
        @info "Ionisation (Lucas-Saelee) Test"
        T = 0.
        meas, swarm, params = mc(LucasSaelee(),E=1.,T=T,N=N,t=3e-6)
        sol = pt(LucasSaelee(),E=1.,T=T,lmax=9,N=5000)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
end
