module nanowire

using NanoWiresJulia
using GradientRobustMultiPhysics
using ExtendableGrids
using GridVisualize
using Printf
using SimplexGridFactory
using Triangulate
using TetGen
using DrWatson
using DataFrames

## start with run_watson() --> result goto data directory
## postprocess data with postprocess(; Plotter = PyPlot) --> images go to plots directory

# configure Watson
@quickactivate "NanoWiresJulia" # <- project name
# set parameters that should be included in filename
watson_accesses = ["geometry", "scenario", "shell_x", "stressor_x", "full_nonlin", "nrefs", "femorder", "femorder_P"] 
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
        "femorder_P" => 1,                      # order of the finite element discretisation (polarisation)
        "upscaling" => 1,                       # upscaling of results (does so many extra nrefs for plots)
        "nrefs" => 0,                           # number of uniform refinements before solve
        "avgc" => 2,                            # lattice number calculation method (average case)
        "polarisation" => true,                 # also solve for polarisation
        "fully_coupled" => false,               # parameter for later when we have the full model
        "postprocess" => true,                  # angle calculation, vtk files, cuts
    )
    return params
end

function set_params!(d,kwargs)
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
function main(d = nothing; verbosity = 0, Plotter = nothing, force::Bool = false, kwargs...)

    ## load parameter set
    if d === nothing
        d = get_defaults()
        set_params!(d, kwargs)
    end

    println("***Solving nanowire problem***")

    # TODO: a grid generator in julia would be nice
    # println("***Creating nanowire mesh***")
    # current_directory = pwd()
    # println(current_directory)
    # cd("..")
    # run(`python3 nanowire.py`)
    # cd(current_directory)

    ## set log level
    set_verbosity(verbosity)

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
    gridfile = "grids/nanowire-grid($(geometry[1]),$(geometry[2]),$(geometry[3]),$(geometry[4])).sg" # default: gridfile = "nanowire-grid(30,20,15,2000).sg"
    xgrid = simplexgrid(gridfile)
    xgrid = uniform_refine(xgrid,nrefs)
    @show xgrid


    ###########################
    ### PROBLEM DESCRIPTION ###
    ###########################

    ## create PDEDescription and add displacement unknown
    Problem = PDEDescription("nanowire bending")
    add_unknown!(Problem; unknown_name = "u", equation_name = "displacement equation")
    subiterations = [[1]]

    ## add Dirichlet boundary data on front
    add_boundarydata!(Problem, 1, [1], HomogeneousDirichletBoundary)

    ## add (nonlinear) operators for displacement equation
    for r = 1 : nregions
        if fully_coupled
            # todo
        else
            add_operator!(Problem, 1, get_displacement_operator(MD.TensorC[r], strainm, eps0[r][1], a[r]; dim = 3, emb = parameters, regions = [r], quadorder = quadorder_D)) 
        end
    end

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
    FEType_D = femorder == 1 ? H1P1{3} : H1P2{3,3} 
    FEType_P = femorder_P == 1 ? H1P1{1} : H1P2{1,3}

    ## call solver
    @unpack nsteps, tres, maxits = d
    Solution, residual = solve_by_embedding(Problem, xgrid, parameters;
                    subiterations = subiterations,
                    nsteps = [nsteps, 1],
                    FETypes = [FEType_D, FEType_P],
                    target_residual = [tres, tres],
                    maxiterations = [maxits, 1])

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
        postprocess(filename; Plotter = Plotter)
    end

    return resultd
end


function get_lattice_misfit_nanowire(lcavg_case, MD::MaterialData, geometry)

    nregions = length(MD.data)
    lc_avg::Array{Float64,1} = zeros(Float64,3)
    eps0::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, nregions)
    r::Array{Int64,1} = geometry[1:3]
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
    set_params!(d, kwargs)

    # load dict from file
    filename = savename(d, "jld2"; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    d = wload(datadir(watson_datasubdir, filename))
    return d
end

function postprocess(filename; Plotter = nothing, cut_levels = "auto", cut_npoints = 100, vol_cut = "auto")

    if typof(filename) <: Dict
        d = filename
    else
        d = wload(datadir(watson_datasubdir, filename))
    end

    ## compute statistics
    @unpack solution, geometry = d
    angle, curvature, dist_bend, farthest_point = compute_statistics(solution[1].FES.xgrid, solution[1], [0,0,geometry[4]])
    d["angle"] = angle
    d["curvature"] = curvature

    # export vtk files
    @unpack polarisation, strainm, upscaling = d
    filename_vtk = savename(d, ""; allowedtypes = watson_allowedtypes, accesses = watson_accesses)
    if polarisation
        writeVTK(datadir(watson_datasubdir, filename_vtk), solution[1], solution[2]; upscaling = upscaling, strain_model = strainm, eps_gfind = 1e-10)
    else
        writeVTK(datadir(watson_datasubdir, filename_vtk), solution[1]; upscaling = upscaling, strain_model = strainm, eps_gfind = 1e-10)
    end

    ## save again
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
    perform_plane_cuts(datadir(watson_datasubdir, filename_cuts), solution, geometry, cut_levels; strain_model = strainm, cut_npoints = cut_npoints, vol_cut = vol_cut, Plotter = Plotter)
end


function perform_plane_cuts(target_folder_cut, Solution, geometry, cut_levels; strain_model = NonlinearStrain3D, cut_npoints = 100, vol_cut = 16, eps_gfind = 1e-11, Plotter = nothing)

    xgrid = Solution[1].FES.xgrid
    diam = geometry[1] + geometry[2]


    ## find three points on the plane z = cut_level and evaluate displacement at points of plane
    @info "Calculating coefficients of plane equations for cuts at levels $(cut_levels)"
    xref = [zeros(Float64,3),zeros(Float64,3),zeros(Float64,3)]
    cells = zeros(Int,3)
    PE = PointEvaluator(Solution[1], Identity)
    CF = CellFinder(xgrid)

    plane_equation_coeffs = Array{Array{Float64,1},1}(undef, length(cut_levels))
    for l = 1 : length(cut_levels)

        ## define plane equation coefficients
        ## 1:3 = normal vector
        ##   4 = - normal vector ⋅ point on plane
        ## find normal vector of displaced plane defined by the three points x[1], x[2] and x[3] 
        cut_level = cut_levels[l]

        x = [[-0.25*diam,-0.25*diam,cut_level],[0.25*diam,-0.25*diam,cut_level],[-0.25*diam,0.25*diam,cut_level]]

        result = deepcopy(x[1])
        for i = 1 : 3
            # find cell
            cells[i] = gFindLocal!(xref[i], CF, x[i]; icellstart = 1, eps = eps_gfind)
            if cells[i] == 0
                cells[i] = gFindBruteForce!(xref[i], CF, x[i])
            end
            @assert cells[i] > 0
            # evaluate displacement
            evaluate!(result,PE,xref[i],cells[i])
            ## displace point
            x[i] .+= result
        end

        plane_equation_coeffs[l] = zeros(Float64,4)
        plane_equation_coeffs[l][1]  = (x[1][2] - x[2][2]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[l][1] -= (x[1][3] - x[2][3]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[l][2]  = (x[1][3] - x[2][3]) * (x[1][1] - x[3][1])
        plane_equation_coeffs[l][2] -= (x[1][1] - x[2][1]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[l][3]  = (x[1][1] - x[2][1]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[l][3] -= (x[1][2] - x[2][2]) * (x[1][1] - x[3][1])
        plane_equation_coeffs[l] ./= sqrt(sum(plane_equation_coeffs[l].^2))
        plane_equation_coeffs[l][4] = -sum(x[1] .* plane_equation_coeffs[l][1:3])
    end

    ## interpolate strain into P1 space
    Strain = FEVector{Float64}("ϵ(u_h)",FESpace{H1P1{6}}(xgrid))
    xCoordinates = xgrid[Coordinates]
    nnodes = size(xCoordinates,2)
    nodevals = nodevalues(Solution[1], Gradient)
    strain = zeros(Float64,6)
    for j = 1 : nnodes
        eval_strain!(strain,view(nodevals,:,j), strain_model)
        for k = 1 : 6
            Strain.entries[(k-1)*nnodes+j] = strain[k]
        end
    end

    ## displace grid
    displace_mesh!(xgrid, Solution[1])

    ## cut displaced grid at plane
    for l = 1 : length(cut_levels)
        cut_level = cut_levels[l]
        @info "Cutting domain at wire-length s = $cut_level with plane equation coefficients $(plane_equation_coeffs[l])"
        @time cut_grid, xgrid_uni, xtrafo! = get_cutgrids(xgrid, plane_equation_coeffs[l]; npoints = cut_npoints, vol_cut = vol_cut)

        # plot boundary-conforming Delaunay cut mesh (suitable for FV)
        @info "Plotting Delaunay cut mesh..."
        gridplot(cut_grid, Plotter = Plotter, title = "Delaunay mesh of cut", fignumber = 1)

        ## interpolate data on uniform cut_grid
        @info "Interpolating data on uniform cut mesh..."
        FES2D = FESpace{H1P1{3}}(xgrid_uni)
        FES2D_ϵ = FESpace{H1P1{6}}(xgrid_uni)
        FES2D_P = FESpace{H1P1{1}}(xgrid_uni)
        CutSolution_u = FEVector{Float64}("u (on 2D cut at z = $(cut_level))", FES2D)
        CutSolution_ϵu = FEVector{Float64}("ϵ(u) (on 2D cut at z = $(cut_level))", FES2D_ϵ)
        CutSolution_P = FEVector{Float64}("P (on 2D cut at z = $(cut_level))", FES2D_P)
        @time interpolate!(CutSolution_u[1], Solution[1]; xtrafo = xtrafo!, not_in_domain_value = NaN, only_localsearch = true, eps = eps_gfind)
        @time interpolate!(CutSolution_ϵu[1], Strain[1]; xtrafo = xtrafo!, not_in_domain_value = NaN, only_localsearch = true, eps = eps_gfind)
        if length(Solution) > 1
            @time interpolate!(CutSolution_P[1], Solution[2]; xtrafo = xtrafo!, not_in_domain_value = NaN, only_localsearch = true, eps = eps_gfind)
        end

        ## write data into csv file
        @info "Writing data into csv file..."
        writeVTK!(target_folder_cut * "cut_$(cut_level)_data.vtu", [CutSolution_u[1],CutSolution_ϵu[1],CutSolution_P[1]]; operators = [Identity, Identity, Identity])
        writeCSV!(target_folder_cut * "cut_$(cut_level)_data.txt", [CutSolution_u[1],CutSolution_ϵu[1],CutSolution_P[1]]; operators = [Identity, Identity, Identity], seperator = "\t")

        ## replacing NaN with 1e30 so that min/max calculation works
        replace!(CutSolution_u.entries, NaN=>1e30)
        replace!(CutSolution_ϵu.entries, NaN=>1e30)
        replace!(CutSolution_P.entries, NaN=>1e30)

        ## plot displacement, strain and polarisation on uniform cut grid
        @info "Plotting data on uniform cut grid..."
        uxmin::Float64 = 1e30
        uxmax::Float64 = -1e30
        uymin::Float64 = 1e30
        uymax::Float64 = -1e30
        uzmin::Float64 = 1e30
        uzmax::Float64 = -1e30
        Pmin::Float64 = 1e30
        Pmax::Float64 = -1e30
        ϵmax = -1e30*ones(Float64,6)
        ϵmin = 1e30*ones(Float64,6)
        nnodes_uni = size(xgrid_uni[Coordinates],2)
        for j = 1 : nnodes_uni
            if abs(CutSolution_u.entries[j]) < 1e10
                if length(Solution) > 1
                    Pmin = min(Pmin,CutSolution_P[1][j])
                    Pmax = max(Pmax,CutSolution_P[1][j])
                end
                uxmin = min(uxmin,CutSolution_u[1][j])
                uymin = min(uymin,CutSolution_u[1][nnodes_uni+j])
                uzmin = min(uzmin,CutSolution_u[1][2*nnodes_uni+j])
                uxmax = max(uxmax,CutSolution_u[1][j])
                uymax = max(uymax,CutSolution_u[1][nnodes_uni+j])
                uzmax = max(uzmax,CutSolution_u[1][2*nnodes_uni+j])
                for k = 1 : 6
                    ϵmax[k] = max(ϵmax[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                    ϵmin[k] = min(ϵmin[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                end
            end
        end
        scalarplot(xgrid_uni, view(CutSolution_u.entries,1:nnodes_uni), Plotter = Plotter; flimits = (uxmin,uxmax), title = "ux on cut", fignumber = 1)
        if isdefined(Plotter,:savefig)
            Plotter.savefig(target_folder_cut * "cut_$(cut_level)_ux.png")
        end
        scalarplot(xgrid_uni, view(CutSolution_u.entries,nnodes_uni+1:2*nnodes_uni), Plotter = Plotter; flimits = (uymin,uymax), title = "uy on cut", fignumber = 1)
        if isdefined(Plotter,:savefig)
            Plotter.savefig(target_folder_cut * "cut_$(cut_level)_uy.png")
        end
        scalarplot(xgrid_uni, view(CutSolution_u.entries,2*nnodes_uni+1:3*nnodes_uni), Plotter = Plotter; flimits = (uzmin,uzmax), title = "uz on cut", fignumber = 1)
        if isdefined(Plotter,:savefig)
            Plotter.savefig(target_folder_cut * "cut_$(cut_level)_uz.png")
        end
        if length(Solution) > 1
            scalarplot(xgrid_uni, CutSolution_P.entries, Plotter = Plotter; flimits = (Pmin,Pmax), title = "Polarisation on cut", fignumber = 1)
            if isdefined(Plotter,:savefig)
                Plotter.savefig(target_folder_cut * "cut_$(cut_level)_P.png")
            end
        end
        for k = 1 : 6
            scalarplot(xgrid_uni, view(CutSolution_ϵu.entries,(k-1)*nnodes_uni+1:k*nnodes_uni), Plotter = Plotter; flimits = (ϵmin[k],ϵmax[k]), title = "ϵu[$k] on cut", fignumber = 1)
            if isdefined(Plotter,:savefig)
                Plotter.savefig(target_folder_cut * "cut_$(cut_level)_ϵ$k.png")
            end
        end
    end
end

end