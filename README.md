NanoWiresJulia
==============

## Working examples

Instantiate this project after cloning or updating via
```
$ julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
```

### Scripts

- `scripts/bimetal_thermal_vs_lattice.jl`: 2D or 3D bimetal simulations (nonlinear elasticity + misfit strain from thermal effects (main_thermal) or lattice mismatch (main_lattice))
- `scripts/bimetal_watson.jl`: runs a set of parametrized bimetal simulations via DrWatson + postprocessing (angle vs. lattice mismatch, 3D bimetal cuts)
- `scripts/nanowire.jl`: 3D nanowire simulation (nonlinear elasticity + polarisation + postprocessing)


## Todo

- grid generator for nanowire in Julia (currently user must provide sg-files, one can be found in grids folder)
- fully coupled elasticity + polarisation model (currently only forward coupling)
- parallelize simulations in run_watson