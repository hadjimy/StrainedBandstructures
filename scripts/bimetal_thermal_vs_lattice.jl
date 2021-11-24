module bimetal

using NanoWiresJulia
using GradientRobustMultiPhysics
using ExtendableGrids
using SimplexGridFactory
using Triangulate
using GridVisualize

## this script has three main functions:
## - main() runs the code with a general misfit strain from whatever source
## - main_thermal() computes a thermal misfit and runs main
## - main_lattice() computes a lattice misfit and runs main
##
## (the model and expected curvature is the same in both cases)

## this calculates a misfit strain from lattice constants and averaging
function main_lattice(; lcavg_case::Int = 2,                # lattice average case (see src/misfit_strains.jl)
                        lattice_factor::Float64 = 0.1,      # lattice factor between material A and B (lc = [5,5*(1+lattice_factor)])
                        C11 = [1200,1200],                  # elasticity tensor parameter C11
                        C44 = [500,500],                    # elasticity tensor parameter C44
                        scale = [100,2000],                 # dimensions of bimetal
                        material_border = 0.5,              # share of material A vs. material B
                        femorder = 2,                       # order of the finite element discretisation
                        nrefinements = 0,                   # number of uniform refinements before solve
                        complicated_model::Bool = true,     # use complicated model
                        solver_embedding::Bool = true,      # use embedding (true) or damping (false) solver ?
                        nsteps::Int = 4,                    # number of embedding steps in embedding solver
                        maxiterations::Int = 100,           # max number of iteration in each embedding step
                        target_residual = 1e-11,            # target residual in each embedding step
                        verbosity = 0,              
                        Plotter = nothing)                  # e.g. Plotter = PyPlot

    ## compute Lame' coefficients μ and λ, and Young's modulus E from coefficients of the elasticity tensor coefficients C11, C44
    ## this part mirrors the setup from the nanowire.jl
    μ = 1e-9 * C44
    λ = 1e-9 * (C11 .- 2*C44)
    E = μ .* (3*λ .+ 2*μ) ./ (λ .+ μ)
    ν = λ ./ (2 * (λ .+ μ))

    ## compute lattice misfit strain
    lc = [5, 5 * ( 1 + lattice_factor )]
    misfit_strain, α = get_lattice_mismatch_bimetal(lcavg_case, [scale[1] * material_border, scale[1] * (1 - material_border)], lc)

    ## run main
    main(; E = E, ν = ν, misfit_strain, α = α, scale = scale, material_border = material_border, complicated_model = complicated_model, nsteps = nsteps, maxiterations = maxiterations, solver_embedding = solver_embedding, target_residual = target_residual, femorder = femorder, nrefinements = nrefinements, verbosity = verbosity, Plotter = Plotter)
end

## this calculates a misfit from thermal conductivity properties
function main_thermal(; ν = [0.3,0.3],                      # elastic numbers of material A and B
                        E = [2.1,1.1],                      # elastic moduli of material A and B
                        ΔT = [580,580],                     # temperature of material A and B
                        α = [1.3e-5,2.4e-5],                # thermal conductivity of material A and B
                        scale = [100,2000],                 # dimensions of bimetal
                        material_border = 0.5,              # share of material A vs. material B
                        femorder = 2,                       # order of the finite element discretisation
                        nrefinements = 0,                   # number of uniform refinements before solve
                        complicated_model::Bool = true,     # use complicated model
                        solver_embedding::Bool = true,      # use embedding (true) or damping (false) solver ?
                        nsteps::Int = 4,                    # number of embedding steps in embedding solver
                        maxiterations::Int = 100,           # max number of iteration in each embedding step
                        target_residual = 1e-11,            # target residual in each embedding step
                        verbosity = 0,              
                        Plotter = nothing)                  # e.g. Plotter = PyPlot

    ## calculate thermal misfit strain, average mechanic is not used
    misfit_strain = ΔT .* α
    α = [0, 0]

    ## run main
    main(; E = E, ν = ν, misfit_strain = misfit_strain, α = α, scale = scale, material_border = material_border, complicated_model = complicated_model, nsteps = nsteps, maxiterations = maxiterations, solver_embedding = solver_embedding, target_residual = target_residual, femorder = femorder, nrefinements = nrefinements, verbosity = verbosity, Plotter = Plotter)
end


## this calculates with user-given misfit strain
function main(; ν = [0.3,0.3],                              # elastic numbers of material A and B
                E = [2.1,1.1],                              # elastic moduli of material A and B
                misfit_strain = [580*1.3e-5,580*2.4e-5],    # misfit strain (use other main functions for thermal/lattice misfits)
                α = [1,1],                                  # averages (use other main functions for thermal/lattice averages)
                scale = [100,2000],                         # dimensions of bimetal
                material_border = 0.5,                      # share of material A vs. material B
                femorder = 2,                               # order of the finite element discretisation
                nrefinements = 0,                           # number of uniform refinements before solve
                complicated_model::Bool = true,             # use complicated model
                solver_embedding::Bool = true,              # use embedding (true) or damping (false) solver ?
                nsteps::Int = 4,                            # number of embedding steps in embedding solver
                maxiterations::Int = 100,                   # max number of iteration in each embedding step
                target_residual = 1e-11,                    # target residual in each embedding step
                verbosity = 0,                      
                Plotter = nothing)                          # e.g. Plotter = PyPlot

    ## set log level
    set_verbosity(verbosity)
    
    ## compute Lame' coefficients μ and λ from ν and E
    μ = E ./ (2  .* (1 .+ ν))
    λ = E .* ν ./ ( (1 .- 2*ν) .* (1 .+ ν))

    ## generate bimetal mesh
    dim::Int = length(scale)
    if dim == 3
        xgrid = bimetal_strip3D(; material_border = material_border, scale = scale)
        @assert femorder in 1:2
        FEType = (femorder == 1) ? H1P1{3} : H1P2{3,3}
    else
        xgrid = bimetal_strip2D(; material_border = material_border, scale = scale)
        FEType = H1Pk{2,2,femorder}
    end
    xgrid = uniform_refine(xgrid,nrefinements)

    ## setup model
    strain_model = (dim == 2) ? NonlinearStrain2D : NonlinearStrain3D # strain model
    stress_tensor = [IsotropicElasticityTensor(λ[1], μ[1], dim), IsotropicElasticityTensor(λ[2], μ[2], dim)] # stress tensor
    emb::Array{Float64,1} = [complicated_model ? 1.0 : 0] # array with embedding parameters for complicated model terms

    ## generate problem description
    Problem = PDEDescription("bimetal deformation under misfit strain")
    add_unknown!(Problem; unknown_name = "u", equation_name = "displacement equation")
    add_operator!(Problem, 1, get_displacement_operator(stress_tensor[1], strain_model, misfit_strain[1], α[1]; dim = dim, emb = emb, regions = [1], quadorder = 2*(femorder-1)))
    add_operator!(Problem, 1, get_displacement_operator(stress_tensor[2], strain_model, misfit_strain[2], α[2]; dim = dim, emb = emb, regions = [2], quadorder = 2*(femorder-1)))
    add_boundarydata!(Problem, 1, [1], HomogeneousDirichletBoundary)
    @show Problem

    ## solve system with FEM
    if solver_embedding
        Solution = solve_by_embedding(Problem, xgrid, emb, nsteps = [nsteps], FETypes = [FEType], target_residual = [target_residual], maxiterations = [maxiterations])
    else
        energy = get_energy_integrator(stress_tensor, strain_model, α; dim = dim)
        Solution = solve_by_damping(Problem, xgrid, energy; FETypes = [FEType], target_residual = target_residual, maxiterations = maxiterations)
    end

    ## compute bending statistics (todo: check curvature formula in 3D)
    analytic_curvature = abs(24 * (misfit_strain[2] - misfit_strain[1]) / (scale[1] * (E[1]/E[2] + E[2]/E[1] + 14)))
    curvature = compute_statistics(xgrid, Solution[1], scale)
    @show analytic_curvature, curvature

    ## displace mesh and plot
    p = GridVisualizer(; Plotter = Plotter, layout = (2,1), clear = true, resolution = (800,600))
    scalarplot!(p[1,1], xgrid, view(nodevalues(Solution[1]; abs = true),1,:), levels = 0, colorbarticks = 7, xlimits = [0, scale[2]+10], ylimits = [-100,scale[1]])
    vectorplot!(p[1,1], xgrid, evaluate(PointEvaluator(Solution[1], Identity)), spacing = [scale[2]/10,scale[1]/2], clear = false, title = "u_h (abs + quiver)")
    xgrid_displaced = displace_mesh(xgrid, Solution[1])
    gridplot!(p[2,1], xgrid_displaced, linewidth = "1", title = "displaced mesh")
end

function get_lattice_mismatch_bimetal(lcavg_case, geometry, lc)
    r::Array{Float64,1} = geometry[1:2]
    a::Array{Float64,1} = zeros(Float64,2)

    A_core = r[1]
    A_stressor = r[2]
    lc_avg = (lc[1]*A_core + lc[2]*A_stressor)/(A_core + A_stressor)

    for j = 1 : 2
        if lcavg_case == 1
            a[j] = (lc[j] - lc[1])/lc[1]
        elseif lcavg_case == 2
            a[j] = (lc[j] - lc_avg)/lc[j]
        elseif lcavg_case == 3
            a[j] = (lc[j] - lc_avg)/lc_avg
        end
    end

    return a .* (1 .+ a./2), a
end

end