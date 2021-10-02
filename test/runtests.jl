using MCSim, FVB
using Test

N = 10000
@testset "MCSim.jl" begin
    @testset "Maxwell" begin
        meas, swarm = mc(MaxwellModel(),E=1.,T=0.,N=N)
        sol = pt(MaxwellModel(),E=1.,T=0.)
        @test sol.ε̄- sol.ε̄*0.01 <= swarm.bε <= sol.ε̄+ sol.ε̄*0.01 
        @test sol.DT - sol.DT *0.01 <= swarm.bDt <= sol.DT + sol.DT *0.01 
        @test sol.DT - sol.DT *0.01 <= swarm.bDt <= sol.DT + sol.DT *0.01 
    end
    @testset "ReidRamp" begin
        meas, swarm = mc(ReidRamp(),E=1.,T=0.,N=N)
        sol = pt(ReidRamp(),E=1.,T=0.)
        @test sol.ε̄- sol.ε̄*0.01 <= swarm.bε <= sol.ε̄+ sol.ε̄*0.01 
        @test sol.DT - sol.DT *0.01 <= swarm.bDt <= sol.DT + sol.DT *0.01 
        @test sol.DT - sol.DT *0.01 <= swarm.bDt <= sol.DT + sol.DT *0.01 
    end
end
