using MCSim

cell = CELL(Δε_detector=0.9, thresh=0.1, z=5e-9, x=1, y=1,r=1e-11)

model = MCMODEL(n=1e5,type=cell, v_0 = GaussianDirectionVelocity(ε=1e2), r_0 = InitSpacial(0,[0,0,1]), n0 = UniformDensity(1e28), En0 = ZeroField(), ion_type = DiscreteIonisation(1.), meas=[:SP_z],optimise_spatial=false)

m, s, p = mc(model,Ar,verbose=false);

diff_SP = zeros(length(p.zbins)-1)
for i in eachindex(diff_SP)
    diff_SP[i] = (m.measg.SP_εs[i+1]-m.measg.SP_εs[i])/(p.zbins[i+1]-p.zbins[i])
end

model2 = MCMODEL(n=1e5,type=cell, v_0 = GaussianDirectionVelocity(ε=120), r_0 = InitSpacial(0,[0,0,1]), n0 = UniformDensity(1e28), En0 = ZeroField(), ion_type = DiscreteIonisation(1.), meas=[:SP_z],optimise_spatial=false)

m2, s2, p2 = mc(model2,Ar,verbose=false);

diff_SP2 = zeros(length(p.zbins)-1)
for i in eachindex(diff_SP)
    diff_SP2[i] = (m2.measg.SP_εs[i+1]-m2.measg.SP_εs[i])/(p2.zbins[i+1]-p2.zbins[i])
end
