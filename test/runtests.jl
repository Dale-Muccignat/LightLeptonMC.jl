using Distributed
addprocs()
using MCSim, FVB 
using Test


N = 100000
tol = 0.01
function withinTol(y,ŷ,tol=0.01)
    @test y - y*tol <= ŷ<= y+y*tol
end
@testset "MCSim.jl" begin
    @testset "Maxwell" begin
        @info "Elastic (Maxwell) Test"
        meas, swarm = mc(MaxwellModel(),E=1.,T=0.,N=N,lmax=9)
        sol = pt(MaxwellModel(),E=1.,T=0.)
        withinTol(sol.ε̄,swarm.bε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "ReidRamp" begin
        @info "Inelastic (ReidRamp) Test"
        meas, swarm = mc(ReidRamp(),E=1.,T=0.,N=N,lmax=9)
        sol = pt(ReidRamp(),E=1.,T=0.)
        withinTol(sol.ε̄,swarm.bε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "ReidRamp - Temperature" begin
        @info "Temperature (ReidRamp) Test"
        T = 298.
        meas, swarm = mc(ReidRamp(),E=1.,T=T,N=N,lmax=9)
        sol = pt(ReidRamp(),E=1.,T=T)
        withinTol(sol.ε̄,swarm.bε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "Lucas-Saelee" begin
        @info "Ionisation (Lucas-Saelee) Test"
        T = 0.
        meas, swarm = mc(LucasSaelee(),E=1.,T=T,N=N,lmax=9)
        sol = pt(LucasSaelee(),E=1.,T=T)
        withinTol(sol.ε̄,swarm.bε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
end
