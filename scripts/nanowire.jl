module nanowire

using NanoWiresJulia
using GradientRobustMultiPhysics
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using Printf
using SimplexGridFactory
using Triangulate
using TetGen
using DrWatson
using DataFrames
using Pardiso

## start with run_watson() --> result goto data directory
## postprocess data with postprocess(; Plotter = PyPlot) --> images go to plots directory

# configure Watson
@quickactivate "NanoWiresJulia" # <- project name
# set parameters that should be included in filename
#watson_accesses = ["mstruct", "geometry", "scenario", "shell_x", "stressor_x", "full_nonlin", "nrefs", "femorder", "femorder_P"]
watson_accesses = ["geometry", "nrefs", "femorder", "stressor_x", "uniform_grid", "interface_refinement", "z_levels_dist"]
watson_allowedtypes = (Real, String, Symbol, Array, DataType)
watson_datasubdir = "nanowire"


## get default parameters
function get_defaults()
    params = Dict(
        "shell_x" => 0.3,                       # x value for x-dependent shell material
        "stressor_x" => 0.5,                    # x value for x-dependent stressor material
        "strainm" => NonlinearStrain3D,         # strain model
        "full_nonlin" => true,                  # use complicated model (ignored if linear strain is used)
        "use_emb" => true,                      # use embedding (true) or damping (false) solver ?
        "nsteps" => 4,                          # number of embedding steps in embedding solver
        "maxits" => 10,                         # max number of iteration in each embedding step
        "tres" => 1e-12,                        # target residual in each embedding step
        "geometry" => [30, 20, 15, 2000],       # dimensions of nanowire
        "scenario" => 1,                        # scenario number that fixes materials for core/shell/stressor
        "mb" => 0.5,                            # share of material A vs. material B
        "mstruct" => ZincBlende001,             # material structure type
        "femorder" => 1,                        # order of the finite element discretisation (displacement)
        "use_lowlevel_solver" => true,          # use new implementation based on low level structures (should be faster)
        "femorder_P" => 1,                      # order of the finite element discretisation (polarisation)
        "nrefs" => 0,                           # number of uniform refinements before solve
        "avgc" => 2,                            # lattice number calculation method (average case)
        "polarisation" => true,                 # also solve for polarisation
        "fully_coupled" => false,               # parameter for later when we have the full model
        "postprocess" => false,                 # angle calculation, vtk files, cuts
        "linsolver" => ExtendableSparse.MKLPardisoLU, # linear solver (try ExtendableSparse.MKLPardisoLU or ExtendableSparse.LUFactorization)
        "uniform_grid" => true,                 # uniform grid in z direction or nonuniform grid with local refinement at cut levels
        "interface_refinement" => false,        # enables a finer grid at material interface for each cross-section
        "z_levels_dist" => 100,                 # distance between z-levels is z_levels_dist * 2^(-nrefs)
        "cut_levels_pos" => [0.5],              # position of cut-levels w.r.t. nanowire length. i.e., cut_levels = cut_levels_pos * geometry[4]
        "grid_version" => 1                     # geometry of stressor; either version 1 or 2
    )
    return params
end

function set_params!(d; kwargs...)
    for (k,v) in kwargs
        d[String(k)]=v
    end
    return nothing
end


function get_scenario(scenario, shell_x, stressor_x, materialstructuretype)
    if scenario == 1
        materials = [GaAs,GaAs,AlInAs{stressor_x}] # scenario 1
    elseif scenario == 2
        materials = [GaAs,AlInAs{shell_x},AlInAs{stressor_x}] # scenario 2
    else
        @error "scenario not defined"
    end
    MD = set_data(materials, materialstructuretype)
    return MD
end

# call this function for the full nanowire scenario
# default parameters (see above) can be modified by kwargs
# Plotter is used for postprocessing
function main(d = nothing; verbosity = 0, Plotter = nothing, force::Bool = false, old_grid = false, kwargs...)

    ## load parameter set
    if d === nothing
        d = get_defaults()
        set_params!(d; kwargs...)
    end

    println("***Solving nanowire problem***")

    ## set log level
    #set_verbosity(verbosity)

    ################
    ### SCENARIO ###
    ################
    ## choose material data DofMapTypes
    @unpack shell_x, stressor_x, scenario, mstruct = d
    MD = get_scenario(scenario, shell_x, stressor_x, mstruct)
    nregions = length(MD.data)
    println("SCENARIO DATA")
    println("=============")
    for j = 1 : length(MD.data)
        println("        material[$j] = $(String(get_materialtype(MD.data[j])))")
    end
    println("          structure = $mstruct")

    #######################
    ### MODEL PARAMETER ###
    #######################
    @unpack geometry, fully_coupled, full_nonlin, strainm, avgc, femorder, femorder_P, polarisation = d
    @assert fully_coupled == false "fully coupled model not yet implemented"
    ## setup parameter Array (which is also returned and can be used later to embedd)
    parameters::Array{Float64,1} = [full_nonlin ? 1 : 0]
    quadorder_D = (femorder > 1) ? 2 : 0 # dramatically increases runtime! (3*(femorder-1) would be exact)
    ## compute lattice misfit strain
    eps0, a = get_lattice_misfit_nanowire(avgc, MD, geometry)

    println("        strain type = $(strainm)")
    println("   elasticity model = $(full_nonlin ? "full" : "simplified")")
    println("      misfit strain = $(eps0)")
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    println("           filename = $filename\n")

    # check if result files already exist
    if isfile(datadir(watson_datasubdir, filename)) && !force
        @warn "Result files already exist, rerun with force = true if needed"
        return nothing
    end

    ############
    ### GRID ###
    ############
    @unpack nrefs = d
    ## load mesh, expecting the following regions
    ##              regions [1,2,3] = [core, shell, stressor]
    ##      boundaryregions [1,2,3] = [core front, shell boundary, stressor boundary]
    if old_grid
        gridfile = "grids/nanowire-grid($(geometry[1]),$(geometry[2]),$(geometry[3]),$(geometry[4])).sg" # default: gridfile = "nanowire-grid(30,20,15,2000).sg"
        xgrid = simplexgrid(gridfile)
        xgrid = uniform_refine(xgrid,nrefs)
    else
        geometry = Array{Float64,1}(geometry)
        geometry[1] /= sqrt(3)  # core hexagon side
        geometry[2] /= sqrt(3)  # shell hexagon side
                                # geometry[3] is the stressor width
        if d["uniform_grid"] == true
            d["cut_levels"] = nothing
        else
            d["cut_levels"] = d["cut_levels_pos"]*d["geometry"][4]   # z-level locations of cuts

        end
        if d["interface_refinement"] == true
            α = geometry[3]/2
        else
            α = nothing
        end
        xgrid = nanowire_tensorgrid(; scale = geometry, nrefs = nrefs, cut_levels = d["cut_levels"], α = α,
            z_levels_dist = d["z_levels_dist"], version = d["grid_version"],
            corner_refinement = false, manual_refinement = true)
        #xgrid = nanowire_tensorgrid(; scale = geometry, nrefs = nrefs, cut_levels = d["cut_levels"], α = α, z_levels_dist = d["z_levels_dist"], version = d["grid_version"])
    end
    #xgrid = nanowire_grid(; scale = geometry)
    #gridplot(xgrid; Plotter=Plotter)
    @show xgrid


    ###########################
    ### PROBLEM DESCRIPTION ###
    ###########################

    ## create PDEDescription and add displacement unknown
    Problem = PDEDescription("nanowire bending")
    add_unknown!(Problem; unknown_name = "u", equation_name = "displacement equation")
    subiterations = [[1]]

    ## add Dirichlet boundary data on front
    if d["interface_refinement"] == true
        regions_bc = [6,7,8] # 6 = core bottom, 7 = shell bottom, 8 = stressor bottom
    else
        regions_bc = [4,5,6]  # 4 = core bottom, 5 = shell bottom, 6 = stressor bottom
    end
    add_boundarydata!(Problem, 1, regions_bc, HomogeneousDirichletBoundary)

    ## add (nonlinear) operators for displacement equation
    for r = 1 : nregions
        if fully_coupled
            # todo
        else
           # add_operator!(Problem, 1, get_displacement_operator(MD.TensorC[r], strainm, eps0[r][1], a[r]; dim = 3, emb = parameters, regions = [r], bonus_quadorder = quadorder_D))
        end
    end
    add_operator!(Problem, 1, get_displacement_operator_new(MD.TensorC, strainm, eps0, a; dim = 3, emb = parameters, regions = 1:nregions, bonus_quadorder = quadorder_D))

    ## add (linear) operators for polarisation equation
    if polarisation
        add_unknown!(Problem; unknown_name = "V_P", equation_name = "polarisation potential equation")
        k0 = 8.854e-3 # 8.854e-12 C/(V m) = 1e9*8.854e-12 C/(V nm)
        kr::Array{Float64,1} = [MD.data[1].kappar, MD.data[2].kappar, MD.data[end].kappar]
        subiterations = fully_coupled ? [[1,2]] : [[1], [2]]

        ## add nonlinear operators for displacement
        for r = 1 : nregions
            if fully_coupled
                # todo
            else
                add_operator!(Problem, [2,1], get_polarisation_from_strain_operator(MD.TensorE[r], strainm, eps0[r][1]; dim = 3, regions = [r]))
                add_operator!(Problem, [2,2], get_polarisation_laplacian_operator(; κ = k0 * kr[r], dim = 3, regions = [r]))
            end
        end
    end
    @show Problem

    ##############
    ### SOLVER ###
    ##############

    ## choose finite element type for displacement (vector-valued) and polarisation (scalar-valued)
    if femorder == 1
        FEType_D = H1P1{3}
        FEType_P = H1P1{1}
    elseif femorder == 2
        FEType_D = H1P2{3,3}
        FEType_P = H1P2{1,3}
    elseif femorder == 3
        FEType_D = H1P3{3,3}
        FEType_P = H1P3{1,3}
    end

    ## call solver
    @unpack nsteps, tres, maxits, linsolver, use_lowlevel_solver = d
    if (use_lowlevel_solver)
        DisplacementOperator = PDEDisplacementOperator(MD.TensorC, strainm, eps0, a, parameters, 3) #get_displacement_operator_new(MD.TensorC, strainm, eps0, a; dim = 3, emb = parameters, regions = 1:nregions, bonus_quadorder = quadorder_D)
        PolarisationOperator = PDEPolarisationOperator(MD.TensorE, strainm, eps0, k0 * kr, 3)
        Solution, residual = solve_lowlevel(xgrid,
                                Problem.BoundaryOperators,
                                Problem.GlobalConstraints,
                                DisplacementOperator,
                                PolarisationOperator,
                                parameters;
                                linsolver = linsolver,
                                nsteps = [nsteps, 1],
                                FETypes = [FEType_D, FEType_P],
                                target_residual = [tres, tres],
                                solve_polarisation = polarisation,
                                coupled = false,
                                maxiterations = [maxits, 1])
    else
        if polarisation
            Solution, residual = solve_by_embedding(Problem, xgrid, parameters;
                            subiterations = subiterations,
                            nsteps = [nsteps, 1],
                            linsolver = linsolver,
                            FETypes = [FEType_D, FEType_P],
                            target_residual = [tres, tres],
                            maxiterations = [maxits, 1])
        else
            Solution, residual = solve_by_embedding(Problem, xgrid, parameters;
                            subiterations = subiterations,
                            nsteps = [nsteps],
                            linsolver = linsolver,
                            FETypes = [FEType_D],
                            target_residual = [tres],
                            maxiterations = [maxits])
        end
    end

    ## add computed data to dict
    resultd = deepcopy(d)
    resultd["eps0"] = eps0
    resultd["a"] = a
    resultd["solution"] = Solution
    resultd["residual"] = residual
    wsave(datadir(watson_datasubdir, filename), resultd)

    ###################
    ### POSTPROCESS ###
    ###################
    if d["postprocess"]
        d = postprocess(filename; Plotter = Plotter, cut_levels = geometry[4]/2, simple_cuts = true)
        resultd["angle"] = d["angle"]
        resultd["curvature"] = d["curvature"]
        wsave(datadir(watson_datasubdir, filename), resultd)
    else
        @info "skipping postprocessing... start it manually with: nanowire.postprocess(\"$filename\"; Plotter = PyPlot)"
    end

    return resultd
end


function get_lattice_misfit_nanowire(lcavg_case, MD::MaterialData, geometry)

    nregions = length(MD.data)
    lc_avg::Array{Float64,1} = zeros(Float64,3)
    eps0::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, nregions)
    r::Array{Float64,1} = geometry[1:3]
    a::Array{Float64,1} = zeros(Float64,3)

    lc = Array{Array{Float64,1}}(undef, nregions)
    for j = 1 : nregions
        lc[j] = MD.data[j].LatticeConstants
    end

    if nregions == 2
        # equalibrium strain for nanowire
        A_core = r[1]*r[2]
        A_stressor = r[1]*r[3]
        for j = 1 : 3
            #lc_avg[j] = (lc[1][j]*A_core + lc[2][j]*A_stressor)/(A_core + A_stressor)
            lc_avg[j] = 2*lc[1][j]*lc[2][j]/(lc[1][j] + lc[2][j])
        end
    elseif nregions == 3
        # equilibrium strain for nanowire
        A_core = 3*sqrt(3)/2 * r[1]^2
        A_shell = 3*sqrt(3)/2 * (r[2]^2 + 2*r[1]*r[2])
        A_stressor = sqrt(3)/2 * r[3] * (7*(r[1]+r[2]) + 3*r[3])
        for j = 1 : 3
            lc_avg[j] = (lc[1][j]*A_core + lc[2][j]*A_shell + lc[3][j]*A_stressor)/(A_core + A_shell + A_stressor)
        end
    end

    for region = 1 : nregions
        eps0[region] = zeros(Float64,3)
        for j = 1 : 3
            if lcavg_case == 1
                a[j] = (lc[region][j] - lc[1][j])/lc[1][j]
            elseif lcavg_case == 2
                a[j] = (lc[region][j] - lc_avg[j])/lc[region][j]
            elseif lcavg_case == 3
                a[j] = (lc[region][j] - lc_avg[j])/lc_avg[j]
            end
            eps0[region][j] = a[j] * (1 + a[j]/2)
        end
    end

    return eps0, a
end

## load a result dictionary (to browse results in julia)
function load_data(d = nothing; kwargs...)
    ## complete parameter set
    if d === nothing
        d = get_defaults()
    end
    set_params!(d; kwargs...)

    # load dict from file
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    d = wload(datadir(watson_datasubdir, filename))
    return d
end

function export_vtk(d = nothing; upscaling = 0, kwargs...)
    d = load_data(d; kwargs...)
    filename_vtk = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    solution = d["solution"]
    repair_grid!(solution[1].FES.xgrid)
    exportVTK(datadir(watson_datasubdir, filename_vtk), solution[1]; P0strain = true, upscaling = upscaling, strain_model = d["strainm"], eps0 = d["eps0"])
end

function postprocess(filename = nothing; Plotter = nothing, export_vtk = true, cut_levels = "auto", simple_cuts = true, cut_npoints = 200, vol_cut = "auto", eps_gfind = 1e-12, upscaling = 0, kwargs...)

    if typeof(filename) <: Dict
        d = filename
        filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    elseif filename === nothing
        d = load_data(; kwargs...)
        filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    else
        d = wload(datadir(watson_datasubdir, filename))
    end

    ## compute statistics
    @unpack solution, geometry = d
    repair_grid!(solution[1].FES.xgrid)
    angle, curvature, dist_bend, farthest_point = compute_statistics(solution[1].FES.xgrid, solution[1], geometry, eltype(solution[1].FES))
    d["angle"] = angle
    d["curvature"] = curvature

    # export vtk files
    if export_vtk == true
        @unpack polarisation, strainm, eps0 = d
        filename_vtk = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
        if polarisation
            exportVTK(datadir(watson_datasubdir, filename_vtk), solution[1], solution[2]; P0strain = true, upscaling = upscaling, strain_model = strainm, eps0 = eps0)
        else
            exportVTK(datadir(watson_datasubdir, filename_vtk), solution[1]; P0strain = true, upscaling = upscaling, strain_model = strainm, eps0 = eps0)
        end
    end

    ## save again
    @show filename
    wsave(datadir(watson_datasubdir, filename), d)

    ## compute cuts (only exported as vtu and png files)
    @unpack strainm, nrefs = d
    if cut_levels == "auto"
        cut_levels = [0.3,0.4,0.5,0.6,0.7] * geometry[4]
    end
    if vol_cut == "auto"
        vol_cut = 4.0^(2-nrefs)
    end
    filename_cuts = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses) * "_CUTS/"
    mkpath(datadir(watson_datasubdir, filename_cuts))
    diam = geometry[1] + geometry[2]
    plane_points = [[-0.25*diam,-0.25*diam],[0.25*diam,-0.25*diam],[-0.25*diam,0.25*diam]] # = [X,Y] coordinates of the three points that define the cut plane
    if simple_cuts # needs grid that triangulates cut_levels
        perform_simple_plane_cuts(datadir(watson_datasubdir, filename_cuts), solution, plane_points, cut_levels; eps0 = d["eps0"], eps_gfind = eps_gfind, cut_npoints = cut_npoints, only_localsearch = true, strain_model = strainm, Plotter = Plotter, upscaling = upscaling)
    else
        #perform_plane_cuts(datadir(watson_datasubdir, filename_cuts), solution, plane_points, cut_levels; strain_model = strainm, cut_npoints = cut_npoints, vol_cut = vol_cut, Plotter = Plotter)
    end

    return d

end

end
