using Test

using MCSim, Gases, LXCatCompute

@testset "Individual Checks" begin
    @testset "Z bin changes" begin
        cell = CELL(z=5mm,z_n=21)
        model = MCMODEL(type=cell,meas=(:tz,),r_0=InitSpatial(1.2 *mm,[0,0,1]),En0=UniformField(10Td),v_0=InitVelocity(1eV,[0,0,1]))
        params = initParams(Gas(Ne),model)
        meas_inner = initMeas(model,params)
        particle = PART(model,params)
        dt = MCSim.MC.getZDt(model.En0,particle,model,params)
        updateParticle!(dt.time,particle,model.En0,model.n0,model,params,meas_inner,false)
        # MCSim.MC.event!(dt,particle,model,params,meas_inner,false)
        @test params.zbins[particle.zbin+1] <
particle.r[3] < params.zbins[particle.zbin+2]

        model = MCMODEL(type=cell,meas=(:tz,),r_0=InitSpatial(1.2 *mm,[0,0,1]),En0=UniformField(150Td),v_0=InitVelocity(1eV,[0,0,-1]))
        params = initParams(Gas(Ne),model)
        meas_inner = initMeas(model,params)
        particle = PART(model,params)
        dt = MCSim.MC.getZDt(model.En0,particle,model,params)
        updateParticle!(dt.time,particle,model.En0,model.n0,model,params,meas_inner,false)
        # MCSim.MC.event!(dt,particle,model,params,meas_inner,false)

        @test params.zbins[particle.zbin+1] < particle.r[3] < params.zbins[particle.zbin+2]


        model = MCMODEL(type=cell,meas=(:tz,),r_0=InitSpatial(1.2 *mm,[0,0,1]),En0=UniformField(10Td),v_0=InitVelocity(100eV,[0,0,-1]))
        params = initParams(Gas(Ne),model)
        meas_inner = initMeas(model,params)
        particle = PART(model,params)
        dt = MCSim.MC.getZDt(model.En0,particle,model,params)
        updateParticle!(dt.time,particle,model.En0,model.n0,model,params,meas_inner,false)
        # MCSim.MC.event!(dt,particle,model,params,meas_inner,false)

        @test params.zbins[particle.zbin-1] < particle.r[3] < params.zbins[particle.zbin]
    end

    @testset "Collision angle" begin
        model = MCMODEL()
        params = initParams(MaxwellModel(),model)
        meas_inner = initMeas(model,params)
        particle = PART(model,params)

        # coschi = 1-2*rand()
        results_angle = []
        results_mag = []
        for _ in 1:100
            coschi = params.splat_gas.dcs[cs_no][particle.cbin](rand())::Float64
            #coschi = -1
            sinchi = sqrt(1-coschi^2)
            phi = rand()*2*pi
            vecs = [rand(3) for _ in 1:100] |> (x) -> x./norm.(x)
            v2 = [cos(phi)*sinchi,sin(phi)*sinchi,coschi]
            new_vecs = [MCSim.MC.shift(v2,vec) for vec in vecs]
            coschis = [dot(new_vec,vec)/norm(vec)/norm(new_vec) for (new_vec,vec) in zip(new_vecs,vecs)]
            push!(results_angle, minimum(coschis .≈ coschi))
            norms = [norm(vec) ≈ norm(new_vec) for (new_vec,vec) in zip(new_vecs,vecs)]
        end
        @test minimum(results_angle)
        @test minimum(results_mag)

    end
end

