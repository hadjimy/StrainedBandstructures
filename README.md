# StrainedBandstructures

## Working examples

Instantiate this project after cloning or updating via

```
$ julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
```

For example to run the simulation for a GaAs/In0.5Al0.5As nanowire with 50 nm thick core (GaAs) and 10 nm wide shell/stressor (In0.5Al0.5As) type:

```
julia> include("scripts/nanowire.jl")
julia> using PyPlot
julia> nanowire.main(geometry=[30,20,10,2000],stressor_x=0.5,Plotter=PyPlot,force=true,postprocess=true)
```

The code by default uses linear FEM elements without any grid refinement and a ZB-(001) crystal orientation. 
These settings (and many others) can be tuned by changing the papameters in `scripts/nanowire.jl`.

The above commands will generate a .jld2 file with the solution data and two .vtu files with the bent and unbent nanowire that can be imported in Paraview for visualizing the solution.
In addition it will create cross-section cuts and output the elastic strain and other solution derivatives on cross-sections perpedicular to the bending direction (by default a middle cross-section is chosen).

### Scripts

- `scripts/nanowire.jl`: 3D nanowire simulation (nonlinear elasticity + polarisation + postprocessing)
- `scripts/bimetal_watson.jl`: runs a set of parametrized bimetal simulations via DrWatson + postprocessing (angle vs. lattice mismatch, 3D bimetal cuts)
- `scripts/bimetal_thermal_vs_lattice.jl`: 2D or 3D bimetal simulations (nonlinear elasticity + misfit strain from thermal effects (main_thermal) or lattice mismatch (main_lattice))


## Todo

- [Done, 16/03/2022] grid generator for nanowire in Julia (currently user must provide sg-files, one can be found in grids folder)
- fully coupled elasticity + polarisation model (currently only forward coupling)
- parallelize simulations in run_watson
- regression tests
