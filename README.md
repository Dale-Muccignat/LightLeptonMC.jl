# MCSim

A flexible Monte Carlo simulation package for swarm and Liquid Micro-Jet applications.
## Install
To install this package, several dependancies are required and are handled through a local registry. First, add the local registry:
```julia
] registry add git@github.com:Dale-Muccignat/localregistry2.git
```
Then add MCSim
```julia
] add MCSim
```
## Swarm

To use, initiate a `Gas` struct from `Gases.jl`:

```julia
gas = MaxwellModel()
```

and then run:
```julia
meas, swarm_meas = mc(gas,E=1.,T=0.)
```

where `E` is the applied electric field in Townsend units and `T` is the gas Temperature in Kelvin units. `swarm_meas` contains common transport coefficients such as:

- bε: mean energy in eV
- bW: bulk drift velocity in m/s
- bDl: bulk longnitude diffusion
- bDt: bulk transverse diffusion

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

