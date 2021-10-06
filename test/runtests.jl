using MCSim
using Test


N = 1000000
tol = 0.01
function withinTol(y,ŷ,tol=0.01)
    @test y - y*tol <= ŷ<= y+y*tol
end
@testset "MCSim.jl" begin
    @testset "Maxwell" begin
        @info "Elastic (Maxwell) Test"
        meas, swarm = mc(MaxwellModel(),E=1.,T=0.,N=N)
        sol = pt(MaxwellModel(),E=1.,T=0.,lmax=9,N=500)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "ReidRamp" begin
        @info "Inelastic (ReidRamp) Test"
        meas, swarm = mc(ReidRamp(),E=1.,T=0.,N=N)
        sol = pt(ReidRamp(),E=1.,T=0.,lmax=9)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "ReidRamp - Temperature" begin
        @info "Temperature (ReidRamp) Test"
        T = 298.
        meas, swarm = mc(ReidRamp(),E=1.,T=T,N=N)
        sol = pt(ReidRamp(),E=1.,T=T,lmax=9)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
    @testset "Lucas-Saelee" begin
        @info "Ionisation (Lucas-Saelee) Test"
        T = 0.
        meas, swarm = mc(LucasSaelee(),E=1.,T=T,N=N)
        sol = pt(LucasSaelee(),E=1.,T=T,lmax=9)
        withinTol(sol.ε̄,swarm.ε)
        withinTol(sol.DT,swarm.bDt)
        withinTol(sol.DL,swarm.bDl)
    end
end
