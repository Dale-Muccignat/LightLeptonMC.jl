# MCSim

A flexible Monte Carlo simulation package for swarm and Liquid Micro-Jet applications.

## Swarm

To use, initiate a Gas struct from Gases.jl:

```julia
gas = MaxwellModel()
```

and then run:
```julia
meas, swarm_meas = mc(::Gas,E=1.,T=0.)
```

where `E` is the applied electric field in Townsend units and `T` is the gas Temperature in Kelvin units.

## Liquid Micro-Jet

In order to utilize a Liquid Micro-Jet simulation we first initilize a model and gas with:

```julia
model = MCMODEL(type=LμJ())
gas = MaxwellModel()
```

Then run:

```julia
meas, swarm = mc(model,gas)
```
