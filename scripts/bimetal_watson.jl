module bimetal

using NanoWiresJulia
using GradientRobustMultiPhysics
using ExtendableGrids
using SimplexGridFactory
using Triangulate
using GridVisualize
using DrWatson
using DataFrames
using Pardiso
using ExtendableSparse

## start with run_watson() --> result goto data directory
## postprocess data with postprocess(; Plotter = PyPlot) --> images go to plots directory

# configure Watson
@quickactivate "NanoWiresJulia" # <- project name
# set parameters that should be included in filename
watson_accesses = ["scale", "latfac", "femorder", "full_nonlin", "nrefs", "strainm", "mb", "grid_type"]
watson_allowedtypes = (Real, String, Symbol, Array, DataType)
watson_datasubdir = "bimetal_watson"


function run_watson(; dim = 2, force::Bool = false, generate_vtk = true)

    allparams = Dict(
        "latfac" => Array{Float64,1}(0:0.02:0.2),  # lattice factor between material A and B (lc = [5,5*(1+latfac)])
        "E" => [[1e-6, 1e-6]],                     # elastic moduli of material A and B
        "ν" => [[0.15, 0.15]],                     # Poisson numbers of material A and B
        "strainm" => [dim == 2 ? NonlinearStrain2D : NonlinearStrain3D],          # strain model
        "full_nonlin" => [true],                   # use complicated model (ignored if linear strain is used)
        "use_emb" => [true],                       # use embedding (true) or damping (false) solver ?
        "nsteps" => [4],                           # number of embedding steps in embedding solver
        "maxits" => [100],                         # max number of iteration in each embedding step
        "tres" => [1e-12],                         # target residual in each embedding step
        "scale" => dim == 2 ? [[50,2000], [100, 2000]] : [[50,50,2000], [100,100,2000]],       # dimensions of bimetal
        "mb" => [0.5],                             # share of material A vs. material B
        "femorder" => [2],                         # order of the finite element discretisation
        "upscaling" => [0],                        # upscaling of results (does so many extra nrefs for plots)
        "nrefs" => [1],                            # number of uniform refinements before solve
        "avgc" => [2],                             # lattice number calculation method (average case)
        "linsolver" => ExtendableSparse.MKLPardisoLU, # linear solver (try ExtendableSparse.MKLPardisoLU or ExtendableSparse.LUFactorization)
        "grid_type" => "default",                  # grid options: default, condensator, condensator_tensorgrid
    )

    dicts = dict_list(allparams)
    @info "Starting bimetal simulations..." dicts

    # run the simulations and save data
    for (i, d) in enumerate(dicts)
        run_single(d; force = force, generate_vtk = generate_vtk)
    end
end


# defaults for use with run_single
function get_defaults()
    params = Dict(
        "latfac" => 0.05,  # lattice factor between material A and B (lc = [5,5*(1+latfac)])
        "E" => [1e-6, 1e-6],                     # elastic moduli of material A and B
        "ν" => [0.15, 0.15],                     # Poisson numbers of material A and B
        "scale" => [50,2000],                    # dimensions of bimetal
        "full_nonlin" => true,                   # use complicated model (ignored if linear strain is used)
        "use_emb" => true,                       # use embedding (true) or damping (false) solver ?
        "nsteps" => 4,                           # number of embedding steps in embedding solver
        "maxits" => 100,                         # max number of iteration in each embedding step
        "tres" => 1e-12,                         # target residual in each embedding step
        "mb" => 0.5,                             # share of material A vs. material B
        "femorder" => 2,                         # order of the finite element discretisation
        "upscaling" => 0,                        # upscaling of results (does so many extra nrefs for plots)
        "nrefs" => 1,                            # number of uniform refinements before solve
        "avgc" => 2,                             # lattice number calculation method (average case)
        "linsolver" => ExtendableSparse.MKLPardisoLU, # linear solver (try ExtendableSparse.MKLPardisoLU or ExtendableSparse.LUFactorization)
        "grid_type" => "default",                # grid options: default, condensator, condensator_tensorgrid
    )
    dim = length(params["scale"])
    params["strainm"] = dim == 2 ? NonlinearStrain2D : NonlinearStrain3D
    return params
end
function set_params!(d; kwargs)
    for (k,v) in kwargs
        @info "setting $((k,v))"
        d[String(k)]=v
    end
    d["strainm"] = length(d["scale"]) == 2 ? NonlinearStrain2D : NonlinearStrain3D
    return nothing
end

function get_data(d = nothing; kwargs)
    ## load parameter set
    if d === nothing
        d = get_defaults()
    end
    set_params!(d; kwargs)
    return d
end


function run_single(d = nothing; force::Bool = false, generate_vtk = true, kwargs...)

    d = get_data(d; kwargs)
    @show d

    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    if isfile(datadir(watson_datasubdir, filename)) && !force
        @info "Skipping dataset $filename... (run with force = true to enforce recomputation)"
        return nothing
    else
        @info "Running dataset $filename..."
    end

    ## compute lattice_mismatch
    fulld = copy(d)
    latfac = d["latfac"]
    mb = d["mb"]
    avgc = d["avgc"]
    scale = d["scale"]
    lc = [5, 5 * ( 1 + latfac )]
    misfit_strain, α = get_lattice_mismatch_bimetal(avgc, [scale[1] * mb, scale[1] * (1 - mb)], lc)
    fulld["misfit_strain"] = misfit_strain
    fulld["α"] = α

    ## run simulation
    solution, residual = main(fulld)

    ## save additional data
    fulld["solution"] = solution
    fulld["residual"] = residual

    ## save data to file
    wsave(datadir(watson_datasubdir, filename), fulld)

    if generate_vtk
        export_vtk(fulld; upscaling = fulld["upscaling"], kwargs...)
    end

    return fulld
end

## this calculates with user-given misfit strain
function main(d::Dict; verbosity = 0)

    ## unpack paramers
    @unpack linsolver, latfac, E, ν, misfit_strain, α, full_nonlin, use_emb, nsteps, maxits, tres, scale, mb, femorder, nrefs, strainm, avgc, grid_type = d
    
    ## set log level
    set_verbosity(verbosity)
    
    ## compute Lame' coefficients μ and λ from ν and E
    μ = E ./ (2  .* (1 .+ ν))
    λ = E .* ν ./ ( (1 .- 2*ν) .* (1 .+ ν))

    ## generate bimetal mesh
    dim::Int = length(scale)
    @assert (femorder in 1:3)
    if dim == 3
        # default dirichlet regions
        dirichlet_regions = [3,4] ### DOUBLE CHECK!
        if grid_type == "default"
            #xgrid = bimetal_strip3D(; material_border = mb, scale = scale)
            #xgrid = uniform_refine(xgrid,nrefs)
            xgrid = bimetal_tensorgrid(; scale = scale, nrefs = nrefs)
        elseif grid_type == "condensator"
            xgrid = condensator3D(; scale = scale, d = 10, nrefs = nrefs)
            dirichlet_regions = [1,2,5,6] # core sides and bottoms
        elseif grid_type == "condensator_tensorgrid"
            xgrid = condensator3D_tensorgrid(; scale = scale, d = 10, nrefs = nrefs)
            #dirichlet_regions = [7] # stressor side
            dirichlet_regions = [1,2,3,4,5,6] # core sides and bottoms
        else
            xgrid = bimetal_strip3D_middle_layer(; scale = scale, reflevel = nrefs)
        end

        if femorder == 1
            FEType = H1P1{3}
        elseif femorder == 2
            FEType = H1P2{3,3}
        elseif femorder == 3
            FEType = H1P3{3,3}
        end
    else
        xgrid = bimetal_strip2D(; material_border = mb, scale = scale)
        FEType = H1Pk{2,2,femorder}
        xgrid = uniform_refine(xgrid,nrefs)
        dirichlet_regions = [1]
    end

    ## setup model
    full_nonlin *= strainm <: NonlinearStrain
    emb::Array{Float64,1} = [full_nonlin ? 1.0 : 0] # array with embedding parameters for complicated model terms

    ## generate problem description
    Problem = PDEDescription("bimetal deformation under misfit strain")
    add_unknown!(Problem; unknown_name = "u", equation_name = "displacement equation")
    add_operator!(Problem, 1, get_displacement_operator(IsotropicElasticityTensor(λ[1], μ[1], dim), strainm, misfit_strain[1], α[1]; dim = dim, emb = emb, regions = [1], bonus_quadorder = 2*(femorder-1)))
    add_operator!(Problem, 1, get_displacement_operator(IsotropicElasticityTensor(λ[2], μ[2], dim), strainm, misfit_strain[2], α[2]; dim = dim, emb = emb, regions = [2], bonus_quadorder = 2*(femorder-1)))
    if length(dirichlet_regions) > 0
        add_boundarydata!(Problem, 1, dirichlet_regions, HomogeneousDirichletBoundary)
    else
        add_constraint!(Problem, FixedDofs(1, [1], [0]))
        add_operator!(Problem, [1,1], BilinearForm([NormalFlux, NormalFlux]; factor = 1e10, AT = ON_BFACES, regions = [1]))
    end
    @show Problem

    ## solve system with FEM
    if use_emb
        Solution, residual = solve_by_embedding(Problem, xgrid, emb, nsteps = [nsteps], 
        linsolver = linsolver, FETypes = [FEType], target_residual = [tres], maxiterations = [maxits])
    else
        energy = get_energy_integrator(stress_tensor, strainm, α; dim = dim)
        Solution, residual = solve_by_damping(Problem, xgrid, energy; FETypes = [FEType], linsolver = linsolver, target_residual = tres, maxiterations = maxits)
    end

    return Solution, residual
end

function get_lattice_mismatch_bimetal(avgc, geometry, lc)
    r::Array{Float64,1} = geometry[1:2]
    a::Array{Float64,1} = zeros(Float64,2)

    A_core = r[1]
    A_stressor = r[2]
    lc_avg = (lc[1]*A_core + lc[2]*A_stressor)/(A_core + A_stressor)

    for j = 1 : 2
        if avgc == 1
            a[j] = (lc[j] - lc[1])/lc[1]
        elseif avgc == 2
            a[j] = (lc[j] - lc_avg)/lc[j]
        elseif avgc == 3
            a[j] = (lc[j] - lc_avg)/lc_avg
        end
    end

    return a .* (1 .+ a./2), a
end


## load a result dictionary (to browse results in julia)
function load_data(d = nothing; kwargs)
    ## complete parameter set
    d = get_data(d; kwargs)

    # load dict from file
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    d = wload(datadir(watson_datasubdir, filename))
    return d
end

function export_vtk(d = nothing; upscaling = 0, kwargs...)
    d = load_data(d; kwargs)
    filename_vtk = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    solution = d["solution"]
    repair_grid!(solution[1].FES.xgrid)
    NanoWiresJulia.exportVTK(datadir(watson_datasubdir, filename_vtk), solution[1]; upscaling = upscaling, strain_model = d["strainm"])
end

function export_cuts(; 
    scale = [50, 50, 2000],
    nrefs = 1,
    upscaling = 0,
    eps_gfind = 1e-11,
    mb = 0.5,
    femorder = "max",
    strainm = length(scale) == 3 ? NonlinearStrain3D : NonlinearStrain2D,
    cut_levels = [scale[3] * 0.5],
    Plotter = nothing)

    @info "Exporting cuts..."

    # load all data
    alldata = collect_results(datadir(watson_datasubdir))

    # filter and sort
    df = filter(:scale => ==(scale), alldata)
    df = filter(:nrefs => ==(nrefs), df)
    df = filter(:mb => ==(mb), df)
    df = filter(:strainm => ==(strainm), df)

    if femorder == "max"
        femorder = maximum(df[!,:femorder])
    end
    if upscaling == "auto"
        upscaling::Int = femorder > 1
    end

    df = filter(:femorder => ==(femorder), df)
    @show df

    # compute cuts
    for data in eachrow(df)
        solution = data[:solution]
        strainm = data[:strainm]

        repair_grid!(solution[1].FES.xgrid)
        xgrid_plot = solution[1].FES.xgrid
        
        Solution_plot = nothing
        if upscaling > 0
            xgrid_plot = uniform_refine(xgrid_plot, upscaling; store_parents = true)
            if length(solution) == 1
                FES = FESpace{eltype(solution[1].FES)}(xgrid_plot)
                Solution_plot = FEVector{Float64}("u_h (upscale)", FES)
            else
                FES = [FESpace{eltype(solution[1].FES)}(xgrid_plot), FESpace{eltype(solution[2].FES)}(xgrid_plot)]
                Solution_plot = FEVector{Float64}(["u_h (upscale)", "V_P (upscale)"],[FES[1], FES[2]])
            end
            for j = 1 : length(solution)
                interpolate!(Solution_plot[j], solution[j]; use_cellparents = true, eps = eps_gfind)
            end
        else
            Solution_plot = solution
        end

        filename_cuts = savename(data, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses) * "_CUTS/"
        mkpath(datadir(watson_datasubdir, filename_cuts))
        plane_points = [[0.25*scale[1],0.25*scale[2]],[0.75*scale[1],0.25*scale[2]],[0.25*scale[1],0.75*scale[2]]] # = [X,Y] coordinates of the three points that define the cut plane
        perform_simple_plane_cuts(datadir(watson_datasubdir, filename_cuts), Solution_plot, plane_points, cut_levels; eps_gfind = 1e-10, only_localsearch = true, strain_model = strainm, Plotter = Plotter)
    end
end

function postprocess(; 
    scales = [[50, 2000], [100, 2000]],
    maxlc = [0.1, 0.2],
    nrefs = 1,
    mb = 0.5,
    femorder = "max",
    strainm = length(scales[1]) == 3 ? NonlinearStrain3D : NonlinearStrain2D,
    Plotter = nothing)

    @assert Plotter !== nothing "need a Plotter (e.g. PyPlot)"
    Plotter.close("all")
    @info "Starting postprocessing..."

    # load all data
    alldata = collect_results(datadir(watson_datasubdir))

    # init plot
    fig, (ax1, ax2) = Plotter.subplots(2, 1);
    marker = ["o","s"]
    color = ["red", "blue"]
    legend = []

    for j = 1 : length(scales)

        # filter and sort
        scale = scales[j]
        df = filter(:scale => ==(scale), alldata)
        df = filter(:nrefs => ==(nrefs), df)
        df = filter(:mb => ==(mb), df)
        df = filter(:strainm => ==(strainm), df)
        df = filter(:latfac => <=(maxlc[j]), df)

        if femorder == "max"
            femorder = maximum(df[!,:femorder])
        end
        df = filter(:femorder => ==(femorder), df)
        sort!(df, [:latfac])
        @show df

        # compute curvature etc
        lattice_mismatch = []
        sim_angle = []
        ana_angle = []
        sim_curvature = []
        ana_curvature = []
        for data in eachrow(df)
            solution = data[:solution]
            scale = data[:scale]
        # misfit_strain = data[:misfit_strain]
                ## compute bending statistics (todo: check curvature formula in 3D)
            angle, curvature, dist_bend, farthest_point = compute_statistics(solution[1].FES.xgrid, solution[1], scale)

            # calculate analytic curvature (is it correct for 3D?)
            E = data[:E]
            α = data[:α]
            mb = data[:mb]
            factor = 1/2*(α[2] - α[1])*(2 + α[1] + α[2])
            h1 = data[:scale][1] * mb
            h2 = data[:scale][1] * (1-mb)
            m = h1/h2
            n = E[1]/E[2]
            analytic_curvature = abs( 6.0 * factor * (1+m)^2 / ((h1+h2) * ( 3*(1+m)^2 + (1+m*n)*(m^2+1/(m*n)))) )
            if farthest_point[1] > 0
                analytic_angle = asin(dist_bend/2*analytic_curvature) * 180/π
            else
                analytic_angle = 180 - asin(dist_bend/2*analytic_curvature) * 180/π
            end
            #@info "dist_bend = $(dist_bend)"
            #@info "simulation ===> R = $(1/curvature) | curvature = $curvature | bending angle = $(angle)°"
            #@info "analytic   ===> R = $(1/analytic_curvature) | curvature = $analytic_curvature | bending angle = $(analytic_angle)°"
            push!(lattice_mismatch, data[:latfac])
            push!(sim_angle, angle)
            push!(ana_angle, analytic_angle)
            push!(sim_curvature, curvature)
            push!(ana_curvature, analytic_curvature)
        end
            
        @info "Plotting ..."
        lattice_mismatch *= 100
        ax1.plot(lattice_mismatch, sim_curvature, color = color[j], marker = marker[j])
        ax1.plot(lattice_mismatch, ana_curvature, color = color[j], linestyle = "--")
        ax2.plot(lattice_mismatch, sim_angle, color = color[j], marker = marker[j])
        ax2.plot(lattice_mismatch, ana_angle, color = color[j], linestyle = "--")

        ax1.set_title("Core: $(Int(mb*100))%, Stressor: $(Int((1-mb)*100))%")
        append!(legend, ["simulation, d = $(scale[1])","analytic, d = $(scale[1])"])
        ax1.set_ylabel("curvature")
        ax2.set_ylabel("angle")
        ax2.set_xlabel("lattice mismatch (%)")
    end

    ax1.legend(legend)
    ax2.legend(legend)
    Plotter.savefig("plots/curvature_angle_material_border=$(mb).png")
end

end