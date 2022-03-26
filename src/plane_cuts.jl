

# modification from https://github.com/j-fu/GridVisualize.jl/blob/main/src/common.jl
function tet_x_plane!(ixcoord,ixvalues,iedges,pointlist,node_indices,planeq_values,function_values; tol=0.0)

    # If all nodes lie on one side of the plane, no intersection
    if (mapreduce(a->a< -tol,*,planeq_values) || mapreduce(a->a>tol,*,planeq_values))
        return 0
    end
    # Interpolate coordinates and function_values according to
    # evaluation of the plane equation
    edge_rule = local_celledgenodes(Tetrahedron3D)
    nxs=0
    n1::Int = 0
    n2::Int = 0
    for iedge = 1 : 6
        n1 = edge_rule[1,iedge]
        n2 = edge_rule[2,iedge]
        if planeq_values[n1]*planeq_values[n2]<tol
            nxs+=1
            t= planeq_values[n1]/(planeq_values[n1]-planeq_values[n2])
            for i=1:3
                ixcoord[i,nxs]=pointlist[i,node_indices[n1]]+t*(pointlist[i,node_indices[n2]]-pointlist[i,node_indices[n1]])
            end
            ixvalues[nxs]=function_values[node_indices[n1]]+t*(function_values[node_indices[n2]]-function_values[node_indices[n1]])
            # also remember the local edge numbers that had a succesful cut
            iedges[nxs] = iedge
        end
    end
    return nxs
end

# modification from marching_tetrahedra from https://github.com/j-fu/GridVisualize.jl/blob/main/src/common.jl
function plane_cut(xgrid::ExtendableGrid, plane_equation_coeffs; project_data = [], remesh::Bool = true, vol = 1, tol = 0.0)

    Tv=Float64
    Tp=Float64
    Tf=Int32

    all_ixfaces=Vector{Tf}(undef,0)
    all_ixbedges=Vector{Tf}(undef,0)
    all_ixcoord=Vector{Tp}(undef,0)
    all_ixvalues=Vector{Tv}(undef,0)

    planeq=zeros(4)
    ixcoord=zeros(3,6)
    ixvalues=zeros(6)
    ixbnd=zeros(Bool,6)
    iedges=zeros(Int,6)
    node_indices=zeros(Int32,4)

    plane_equation(plane,coord)=coord[1]*plane[1]+coord[2]*plane[2]+coord[3]*plane[3]+plane[4]

    xCoordinates = xgrid[Coordinates]
    xCellNodes=xgrid[CellNodes]
    nnodes = size(xCoordinates,2)
    func = zeros(Float64,1,nnodes)
    if project_data != []
        for FEB in project_data
            nodevalues!(func, FEB)
        end
    end

    # find nodes along boundary
    edge_on_boundary = zeros(Bool,num_sources(xgrid[GradientRobustMultiPhysics.EdgeNodes]))
    xFaceNodes = xgrid[FaceNodes]
    xFaceCells=xgrid[FaceCells]
    xFaceEdges = xgrid[FaceEdges]
    xCellEdges = xgrid[GradientRobustMultiPhysics.CellEdges]
    for j = 1 : size(xFaceNodes,2)
        if xFaceCells[2,j] == 0
            edge_on_boundary[xFaceEdges[:,j]] .= true
        end
    end

    function pushtris(ns,ixcoord,ixvalues,ixbnd)
        # number of intersection points can be 3 or 4
        if ns>=3
            last_i=length(all_ixvalues)
            for is=1:ns
                # todo: transform points onto z = 0 plane
                @views append!(all_ixcoord,ixcoord[:,is])
                push!(all_ixvalues,ixvalues[is])
            end
            if ixbnd[1] && ixbnd[2]
                append!(all_ixbedges,[last_i+1,last_i+2])
            end
            if ixbnd[2] && ixbnd[3]
                append!(all_ixbedges,[last_i+2,last_i+3])
            end
            if ixbnd[3] && ixbnd[1]
                append!(all_ixbedges,[last_i+3,last_i+1])
            end
            append!(all_ixfaces,[last_i+1,last_i+2,last_i+3])
            if ns==4
                append!(all_ixfaces,[last_i+3,last_i+2,last_i+4])
                if ixbnd[3] && ixbnd[2]
                    append!(all_ixbedges,[last_i+3,last_i+2])
                end
                if ixbnd[2] && ixbnd[4]
                    append!(all_ixbedges,[last_i+2,last_i+4])
                end
                if ixbnd[4] && ixbnd[3]
                    append!(all_ixbedges,[last_i+4,last_i+3])
                end
            end
        end
    end

    start_cell = 0
    for itet=1:size(xCellNodes,2)
        for i=1:4
            node_indices[i]=xCellNodes[i,itet]
        end
        
        @views map!(inode->plane_equation(plane_equation_coeffs,xCoordinates[:,inode]),planeq,node_indices)
        nxs=tet_x_plane!(ixcoord,ixvalues,iedges,xCoordinates,node_indices,planeq,func; tol=tol)
        if nxs >= 3
            start_cell = itet
        end
        # check if cutted edges are boundary edges
        for iedge = 1 : nxs
            ixbnd[iedge] = edge_on_boundary[xCellEdges[iedges[iedge],itet]]
        end
        # save triangles of cut
        pushtris(nxs,ixcoord,ixvalues,ixbnd)
    end

    ## reshape
    xCoordinates = reshape(all_ixcoord,3,Int(length(all_ixcoord)/3))
    nbedges = Int(length(all_ixbedges)/2)
    xBFaceNodes = reshape(all_ixbedges,2,nbedges)

    ## rotate normal around x axis such that n[2] = 0
    ## tan(alpha) = n[2]/n[3]
    alpha = atan(plane_equation_coeffs[2]/plane_equation_coeffs[3])
    #@info "rotating around x-axis with alpha = $alpha"
    plane_equation_coeffs[3] = sin(alpha)*plane_equation_coeffs[2] + cos(alpha)*plane_equation_coeffs[3]
    plane_equation_coeffs[2] = 0

    ## rotate normal around y axis such that n[1] = 0
    ## tan(alpha) = n[2]/n[3]
    beta = -atan(plane_equation_coeffs[1]/plane_equation_coeffs[3])
   # @info "rotating around x-axis with beta = $beta"
    plane_equation_coeffs[3] = -sin(beta)*plane_equation_coeffs[1] + cos(beta)*plane_equation_coeffs[3]
    plane_equation_coeffs[1] = 0

    ## rotate coordinates
    oldcoords = zeros(Float64,3)
    R =  [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)] * [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)]
    for j = 1 : size(xCoordinates,2)
        for k = 1 : 3
            oldcoords[k] = xCoordinates[k,j]
            xCoordinates[k,j] = 0
        end
        for n = 1 : 3, k = 1 : 3
            xCoordinates[n,j] += R[n,k] * oldcoords[k]
        end
    end

    ## restrict coordinates to [x,y] plane
    z = xCoordinates[end,3] # are all the same
    xCoordinates = xCoordinates[1:2,:]

    if remesh
        # problem: boundary nodes are not unique, there are nodes with the same coordinates that we have
        # the following lines look for unique boundary nodes and remaps the BFaceNodes
        bnodes = unique(xBFaceNodes[:])
        node_remap = zeros(Int,length(xCoordinates))
        xCoordinatesB = zeros(Float64,0)
        already_known::Int = 0
        nbnodes::Int = 0
        for node in bnodes
            already_known = 0
            for knode = 1 : nbnodes
                if abs((xCoordinatesB[2*knode-1] - xCoordinates[1,node]).^2 + (xCoordinatesB[2*knode] - xCoordinates[2,node]).^2) < 1e-16
                    already_known = knode
                    break
                end
            end
            if already_known == 0
                append!(xCoordinatesB,xCoordinates[:,node])
                nbnodes += 1
                node_remap[node] = nbnodes
            else
                node_remap[node] = already_known
            end
        end
        @views xBFaceNodes[:] .= node_remap[xBFaceNodes[:]]
        xCoordinatesB = reshape(xCoordinatesB,2,Int(length(xCoordinatesB)/2))
        CenterPoint = sum(xCoordinatesB, dims = 2) ./ size(xCoordinatesB,2)

        # call Grid generator
        xgrid = SimplexGridFactory.simplexgrid(Triangulate;
            points=xCoordinatesB,
            bfaces=xBFaceNodes,
            bfaceregions=ones(Int32,nbedges),
            regionpoints=CenterPoint',
            regionnumbers=[1],
            regionvolumes=[vol])
    else
        ncells = Int(length(all_ixfaces)/3)
        xCellNodes = reshape(all_ixfaces,3,ncells)
        xgrid=ExtendableGrid{Float64,Int32}()
        xgrid[Coordinates] = xCoordinates
        xgrid[CellNodes] =  xCellNodes
        xgrid[CellGeometries] = VectorOfConstants(Triangle2D,ncells);
        xgrid[CellRegions]=ones(Int32,ncells)
        xgrid[BFaceRegions]=ones(Int32,nbedges)
        xgrid[BFaceNodes]=xBFaceNodes
        xgrid[BFaceGeometries]=VectorOfConstants(Edge1D,nbedges)
        xgrid[CoordinateSystem]=Cartesian2D
    end
    
    return xgrid, R, z, start_cell
end



## computes two grids: one boundary conforming Delaunay grid for finite volume methods
## and one uniform non-boundary-conforming one that comprises the cut;
## also a function is returned that transforms 3D coordinates to 2D coordinates on the cut
function get_cutgrids(xgrid, plane_equation_coeffs; npoints = 100, vol_cut = 1.0)
    @time cut_grid, R, z, start_cell = plane_cut(xgrid, plane_equation_coeffs; remesh = true, vol = vol_cut)

    ## define regular grid on cut plane where polarisation should be interpolated
    xmin = minimum(view(cut_grid[Coordinates],1,:))
    xmax = maximum(view(cut_grid[Coordinates],1,:))
    ymin = minimum(view(cut_grid[Coordinates],2,:))
    ymax = maximum(view(cut_grid[Coordinates],2,:))
    @info "Creating uniform grid for bounding box ($xmin,$xmax) x ($ymin,$ymax)"
    h_uni = [xmax - xmin,ymax-ymin] ./ (npoints-1)
    xgrid_uni = simplexgrid(collect(xmin:h_uni[1]:xmax),collect(ymin:h_uni[2]:ymax))

    # mapping from 3D coordinates on cut_grid to 2D coordinates on xgrid_uni
    invR::Matrix{Float64} = inv(R)
    function xtrafo!(x3D,x2D)
        for j = 1 : 3
            x3D[j] = invR[j,1] * x2D[1] + invR[j,2] * x2D[2] + invR[j,3] * z  # invR * [x2D[1],x2D[2],z]
        end
        return nothing
    end

    return cut_grid, xgrid_uni, xtrafo!, start_cell
end


function perform_simple_plane_cuts(target_folder_cut, Solution, plane_points, cut_levels; 
    strain_model = NonlinearStrain3D, 
    export_uniform_data = true,
    cut_npoints = 200, 
    only_localsearch = true,
    eps_gfind = 1e-11,
    Plotter = nothing,
    do_simplecut_plots = Plotter !== nothing,
    do_uniformcut_plots = Plotter !== nothing)

    xgrid = Solution[1].FES.xgrid

    ### first find faces that lie in cut_levels

    # faces with all nodes with the same z-coordinate are assigned to a z-level (otherwise z-level -1 is set)
    xCoordinates = xgrid[Coordinates]
    xCellRegions = xgrid[CellRegions]
    xFaceNodes = xgrid[FaceNodes]
    xFaceCells = xgrid[FaceCells]
    nfaces::Int = num_sources(xFaceNodes)
    z4Faces = zeros(Float64,nfaces)
    nnodes4face = max_num_targets_per_source(xFaceNodes)
    nnodes = size(xCoordinates,2)

    z::Float64 = 0
    for face = 1 : nfaces
        z = xCoordinates[3,xFaceNodes[1,face]]
        for k = 2 : nnodes4face
            if abs(xCoordinates[3,xFaceNodes[k,face]] - z) > eps_gfind
                z = -1
                break;
            end
        end
        z4Faces[face] = z
    end

    ## these are the levels that have faces
    @info "possible z-levels for cut = $(sort(unique(z4Faces))[2:end])"

    ## find three points on the plane z = cut_level and evaluate displacement at points of plane
    xref = [zeros(Float64,3),zeros(Float64,3),zeros(Float64,3)]
    cells = zeros(Int,3)
    PE = PointEvaluator(Solution[1], Identity)
    CF = CellFinder(xgrid)
    plane_equation_coeffs = zeros(Float64,4)

    ## access nodal values of displacement
    nodevals = nodevalues_view(Solution[1])

    ## compute nodal values of gradient of displacement and make a P1 function out of it
    nodevals_gradient = nodevalues(Solution[1], Gradient)
    GlobalInterpolation_∇u = FEVector{Float64}("∇u (global P1 interpolation)", FESpace{H1P1{9}}(xgrid))
    for j = 1 : 9
        view(GlobalInterpolation_∇u.entries,(j-1)*nnodes+1:j*nnodes) .= view(nodevals_gradient,j,:)
    end
    strain = zeros(Float64,6)
    gradient = zeros(Float64,9)

    for l = 1 : length(cut_levels)

        cut_level = cut_levels[l]



        ## find faces for this cut_level
        faces4level = findall(abs.(z4Faces .- cut_level) .< eps_gfind)
        if length(faces4level) < 1
            @warn "found no faces in grid that lie in level = $(cut_level), skipping this cut level..."
            continue
        end
        nodes4level = unique(view(xFaceNodes,:,faces4level))
        node_permute = zeros(Int32,size(xCoordinates,2))
        node_permute[nodes4level] = 1 : length(nodes4level)
        nnodes_cut = length(nodes4level)
        start_cell::Int = xFaceCells[1,faces4level[1]]

        ## define plane equation coefficients and rotation to map 3D coordinates on cut plane to 2D coordinates

        # 1:3 = normal vector
        #   4 = - normal vector ⋅ point on plane
        # find normal vector of displaced plane defined by the three points x[1], x[2] and x[3] 
        x = [[plane_points[1][1],plane_points[1][2],cut_level],[plane_points[2][1],plane_points[2][2],cut_level],[plane_points[3][1],plane_points[3][2],cut_level]]
        result = deepcopy(x[1])
        for i = 1 : 3
            # find cell
            cells[i] = gFindLocal!(xref[i], CF, x[i]; icellstart = start_cell, eps = eps_gfind)
            if cells[i] == 0
                @warn "local search cell search unexpectedly failed, using brute force..."
                cells[i] = gFindBruteForce!(xref[i], CF, x[i])
            end
            @assert cells[i] > 0
            # evaluate displacement
            evaluate!(result,PE,xref[i],cells[i])
            ## displace point
            x[i] .+= result
        end

        plane_equation_coeffs[1]  = (x[1][2] - x[2][2]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[1] -= (x[1][3] - x[2][3]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[2]  = (x[1][3] - x[2][3]) * (x[1][1] - x[3][1])
        plane_equation_coeffs[2] -= (x[1][1] - x[2][1]) * (x[1][3] - x[3][3])
        plane_equation_coeffs[3]  = (x[1][1] - x[2][1]) * (x[1][2] - x[3][2])
        plane_equation_coeffs[3] -= (x[1][2] - x[2][2]) * (x[1][1] - x[3][1])
        plane_equation_coeffs ./= sqrt(sum(plane_equation_coeffs.^2))
        plane_equation_coeffs[4] = -sum(x[1] .* plane_equation_coeffs[1:3])

        ## rotate normal around x axis such that n[2] = 0
        ## tan(alpha) = n[2]/n[3]
        alpha = atan(plane_equation_coeffs[2]/plane_equation_coeffs[3])
        @info "rotating around x-axis with alpha = $alpha"
        plane_equation_coeffs[3] = sin(alpha)*plane_equation_coeffs[2] + cos(alpha)*plane_equation_coeffs[3]
        plane_equation_coeffs[2] = 0

        ## rotate normal around y axis such that n[1] = 0
        ## tan(alpha) = n[2]/n[3]
        beta = -atan(plane_equation_coeffs[1]/plane_equation_coeffs[3])
        @info "rotating around x-axis with beta = $beta"
        plane_equation_coeffs[3] = -sin(beta)*plane_equation_coeffs[1] + cos(beta)*plane_equation_coeffs[3]
        plane_equation_coeffs[1] = 0

        R = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)] * [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)]



        #### SIMPLE CUT GRIDS ###

##        cut_grid = displace_mesh!(xgrid, Solution[1])
        ## map original 3D coordinates onto 2D plane
        xCoordinatesCut3D = zeros(Float64,3,nnodes_cut)
        for j in nodes4level
            ## displaced nodes
            for k = 1 : 3
                xCoordinatesCut3D[k,node_permute[j]] = xCoordinates[k,j] + nodevals[k][j]
            end
        end

        ## construct simple cut grid
        cut_grid=ExtendableGrid{Float64,Int32}()
        cut_grid[Coordinates] = xCoordinatesCut3D
        cut_grid[CellNodes] = node_permute[xFaceNodes[:,faces4level]] # view not possible here
        cut_grid[CellGeometries] = VectorOfConstants{ElementGeometries,Int}(Triangle2D,length(faces4level));
        cut_grid[CellRegions] = xgrid[CellRegions][xgrid[FaceCells][1,faces4level]]
        cut_grid[CoordinateSystem] = Cartesian2D
        cut_grid[BFaceRegions] = ones(Int32,1)
        cut_grid[BFaceNodes] = zeros(Int32,2,0)
        cut_grid[BFaceGeometries] = VectorOfConstants{ElementGeometries,Int}(Edge1D, 0)

        ## interpolate data on cut_grid
        @info "Interpolating data on cut mesh..."
        order_u = get_polynomialorder(eltype(Solution[1].FES), Tetrahedron3D)
        FES2D = FESpace{H1Pk{3,2,order_u}}(cut_grid)
        FES2D_∇u = FESpace{H1P1{9}}(cut_grid)
        FES2D_ϵu = FESpace{H1P1{6}}(cut_grid)
        SimpleCutSolution_u = FEVector{Float64}("u (on 2D cut at z = $(cut_level))", FES2D)
        SimpleCutSolution_∇u = FEVector{Float64}("∇u (on 2D cut at z = $(cut_level))", FES2D_∇u)
        SimpleCutSolution_ϵu = FEVector{Float64}("ϵ(u) (on 2D cut at z = $(cut_level))", FES2D_ϵu)
        
        ## calculate strain from gradient interpolation on cut
        coffset::Int = length(SimpleCutSolution_u.entries)/3
        for j = 1 : nnodes_cut
            for k = 1 : 3
                SimpleCutSolution_u.entries[(k-1)*coffset + j] = nodevals[k][nodes4level[j]]
            end
            for k = 1 : 9 
                SimpleCutSolution_∇u.entries[(k-1)*nnodes_cut + j] = nodevals_gradient[k,nodes4level[j]]
            end
            eval_strain!(strain,view(nodevals_gradient,:,nodes4level[j]), strain_model)
            for k = 1 : 6
                SimpleCutSolution_ϵu.entries[(k-1)*nnodes_cut + j] = strain[k]
            end
        end

        ## 2D coordinates of the simple grid
        cut_grid2D = deepcopy(cut_grid)
        xCoordinatesCutPlane = zeros(Float64,2,nnodes_cut)
        for j = 1 : nnodes_cut
            ## displaced nodes moved to plane
            for n = 1 : 2, k = 1 : 3
                xCoordinatesCutPlane[n,j] += R[n,k] * (xCoordinatesCut3D[k,j])
            end
        end
        cut_grid2D[Coordinates] = xCoordinatesCutPlane
        CF2D = CellFinder(cut_grid2D)

        ## write data into csv file
        @info "Exporting cut data for cut_level = $(cut_level)..."
        kwargs = Dict()
        kwargs[:cellregions] = cut_grid[CellRegions]
        kwargs[:displacement] = view(nodevalues(SimpleCutSolution_u[1], Identity),:,:)
        kwargs[:grad_displacement] = view(nodevalues(SimpleCutSolution_∇u[1], Identity),:,:)
        kwargs[:strain] = view(nodevalues(SimpleCutSolution_ϵu[1], Identity),:,:)
        ExtendableGrids.writeVTK(target_folder_cut * "simple_cut_$(cut_level)_data.vtu", cut_grid; kwargs...)
        
        if do_simplecut_plots
            @info "Plotting data on simple cut grid..."
            nodevals_cut = nodevalues_view(SimpleCutSolution_u[1])
            labels = ["ux","uy","uz"]
            for j = 1 : 3
                scalarplot(cut_grid2D, nodevals_cut[j], Plotter = Plotter; title = "$(labels[j]) on cut", fignumber = 1)
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(target_folder_cut * "simple_cut_$(cut_level)_$(labels[j]).png")
                end
            end
            for k = 1 : 6
                scalarplot(cut_grid2D, view(SimpleCutSolution_ϵu.entries,(k-1)*nnodes_cut+1:k*nnodes_cut), Plotter = Plotter; title = "ϵu[$k] on cut", fignumber = 1)
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(target_folder_cut * "simple_cut_$(cut_level)_ϵ$k.png")
                end
            end
        end


        #### UNIFORM CUT GRIDS ###
        if export_uniform_data


            ## displace grid
            displace_mesh!(xgrid, Solution[1])

            ## map original 3D coordinates onto 2D plane
            xCoordinatesCutPlane = zeros(Float64,3,nnodes_cut)
            for j in nodes4level
                ## displaced nodes moved to plane
                for n = 1 : 3, k = 1 : 3
                    xCoordinatesCutPlane[n,node_permute[j]] += R[n,k] * (xCoordinates[k,j])
                end
            end

            ## restrict coordinates to [x,y] plane
            max_z = maximum(view(xCoordinatesCutPlane,3,:))
            min_z = minimum(view(xCoordinatesCutPlane,3,:))

            ## get z level
            z = (max_z + min_z)/2
            z_error = max_z - min_z
            @info "maximal difference between z coordinates = $(z_error)"
            #invR = inv(R)


            ## define bounding box and uniform cut grid
            xmin = minimum(view(xCoordinatesCutPlane,1,:))
            xmax = maximum(view(xCoordinatesCutPlane,1,:))
            ymin = minimum(view(xCoordinatesCutPlane,2,:))
            ymax = maximum(view(xCoordinatesCutPlane,2,:))
            h_uni = [xmax - xmin,ymax - ymin] ./ (cut_npoints-1)
            Xuni = zeros(Float64,0)
            Yuni = zeros(Float64,0)
            for j = 0 : cut_npoints - 1
                push!(Xuni, h_uni[1]*j)
                push!(Yuni, h_uni[2]*j)
            end
            @info "Creating uniform grid with h_uni = $h_uni for bounding box ($xmin,$xmax) x ($ymin,$ymax)"
            xgrid_uni = simplexgrid(Xuni,Yuni)
            xCoordinatesUni = xgrid_uni[Coordinates]
            nnodes_uni = size(xCoordinatesUni,2)
            
            ## remap simpple grid functions to 2d grid
            FES2D_simple = FESpace{H1Pk{3,2,order_u}}(cut_grid2D)
            FES2D_simple_∇u = FESpace{H1P1{9}}(cut_grid2D)
            FES2D_simple_ϵu = FESpace{H1P1{6}}(cut_grid2D)
            SimpleCut2DSolution_u = FEVector{Float64}("u (on 2D cut at z = $(cut_level))", FES2D_simple)
            SimpleCut2DSolution_∇u = FEVector{Float64}("∇u (on 2D cut at z = $(cut_level))", FES2D_simple_∇u)
            SimpleCut2DSolution_ϵu = FEVector{Float64}("ϵ(u) (on 2D cut at z = $(cut_level))", FES2D_simple_ϵu)
            SimpleCut2DSolution_u.entries .= SimpleCutSolution_u.entries
            SimpleCut2DSolution_∇u.entries .= SimpleCutSolution_∇u.entries
            SimpleCut2DSolution_ϵu.entries .= SimpleCutSolution_ϵu.entries

            ## interpolate data on uniform cut_grid
            @info "Interpolating data on uniform cut mesh..."
            FES2D = FESpace{H1P1{3}}(xgrid_uni)
            FES2D_∇u = FESpace{H1P1{9}}(xgrid_uni)
            FES2D_ϵ = FESpace{H1P1{6}}(xgrid_uni)
            FES2D_P = FESpace{H1P1{1}}(xgrid_uni)
            CutSolution_u = FEVector{Float64}("u (on 2D cut at z = $(cut_level))", FES2D)
            CutSolution_∇u = FEVector{Float64}("∇u (on 2D cut at z = $(cut_level))", FES2D_∇u)
            CutSolution_ϵu = FEVector{Float64}("ϵ(u) (on 2D cut at z = $(cut_level))", FES2D_ϵ)
            CutSolution_P = FEVector{Float64}("P (on 2D cut at z = $(cut_level))", FES2D_P)

            # mapping from 3D coordinates on cut_grid to 2D coordinates on xgrid_uni
            invR::Matrix{Float64} = inv(R)
            function xtrafo!(x3D,x2D)
                for j = 1 : 3
                    x3D[j] = invR[j,1] * (x2D[1] + xmin) + invR[j,2] * (x2D[2] + ymin) + invR[j,3] * z  # invR * [x2D[1],x2D[2],z]
                    ## now x3D is the coordinate of the point in the displaced mesh
                    ## we also need the x3D in the original mesh
                end
                return nothing
            end

            # mapping from 2D coordinates on simple cut_grid to 2D coordinates on xgrid_uni
            function xtrafo2!(x2D_out,x2D_in)
                x2D_out[1] = x2D_in[1] + xmin
                x2D_out[2] = x2D_in[2] + ymin
                return nothing
            end

            interpolate!(CutSolution_u[1], Solution[1]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)

            ## PROBLEM: this interpolates the gradient on the deformed mesh
            ##          how to transform it to the gradient on the original mesh?
            interpolate!(CutSolution_∇u[1], Solution[1]; operator = Gradient, start_cell = start_cell, xtrafo = xtrafo!, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
            if length(Solution) > 1
                interpolate!(CutSolution_P[1], Solution[2]; start_cell = start_cell, xtrafo = xtrafo!, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
            end

            gradmatrix = zeros(Float64,3,3)
            for j = 1 : nnodes_uni
                for k = 1 : 9
                    gradient[k] = CutSolution_∇u.entries[(k-1)*nnodes_uni + j]
                end
                gradmatrix .= (inv([1-gradient[1] gradient[2] gradient[3];
                              gradient[4] 1-gradient[5] gradient[6];
                              gradient[7] gradient[8] 1-gradient[9]]) * [gradient[1] gradient[2] gradient[3];
                              gradient[4] gradient[5] gradient[6];
                              gradient[7] gradient[8] gradient[9]]')
                for k = 1 : 9
                    CutSolution_∇u.entries[(k-1)*nnodes_uni + j] = gradmatrix[k]
                end
                eval_strain!(strain, view(gradmatrix,:), strain_model)
                for k = 1 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes_uni + j] = strain[k]
                end
            end

            ## write material map into txt files
            @info "Writing material map file..."
            filename_material = target_folder_cut * "uniform_cut_$(cut_level)_materialmap.txt"
            io = open(filename_material, "w")
            xCellNodesUni = xgrid_uni[CellNodes]
            nnodes_uni = size(xCoordinatesUni,2)
            xdim = size(xCoordinatesUni,1)
            @printf(io, "WIRE STRESSOR OUTSIDE X Y\n")
            cell::Int = 0
            region::Int = 0
            x3D = zeros(Float64,3)
            xCellRegionsUniform = zeros(Int32,num_cells(xgrid_uni))
            xCellRegionsSimpleCut = cut_grid[CellRegions]

            xtest = zeros(Float64,2)
            for cu = 1 : num_cells(xgrid_uni)
                ## compute center
                fill!(xtest,0)
                for j = 1 : 2, n = 1 : 3
                    xtest[j] += xCoordinatesUni[j, xCellNodesUni[n,cu]]
                end
                xtest ./= 3.0
                xtest[1] += xmin
                xtest[2] += ymin
                cell = gFindLocal!(xref[1], CF2D, xtest; eps = eps_gfind)
                xCellRegionsUniform[cu] = cell == 0 ? 0 : xCellRegionsSimpleCut[cell]
            end
            xgrid_uni[CellRegions] = xCellRegionsUniform

            for n = 1 : nnodes_uni
                for j = 1 : 2
                    xtest[j] = xCoordinatesUni[j,n]
                end
                xtest[1] += xmin
                xtest[2] += ymin
                cell = gFindLocal!(xref[1], CF2D, xtest; eps = eps_gfind)
                region = cell == 0 ? 0 : xCellRegionsSimpleCut[cell]
                if region == 1
                    @printf(io, "1 0 0 ")
                elseif region == 2
                    @printf(io, "0 1 0 ")
                elseif region == 0
                    @printf(io, "0 0 1 ")
                end
                for j = 1 : xdim
                    @printf(io, "%.6f ",xCoordinatesUni[j,n])
                end
                @printf(io, "\n")
            end
            close(io)

            ## write data into vtk file
            @info "Writing data into vtk file..."
            kwargs = Dict()
            kwargs[:cellregions] = xgrid_uni[CellRegions]
            kwargs[:displacement] = view(nodevalues(CutSolution_u[1], Identity),:,:)
            kwargs[:grad_displacement] = view(nodevalues(CutSolution_∇u[1], Identity),:,:)
            kwargs[:strain] = view(nodevalues(CutSolution_ϵu[1], Identity),:,:)
            if length(Solution) > 1
                kwargs[:polarisation] = view(nodevalues(CutSolution_P[1], Identity),:,:)
            end
            ExtendableGrids.writeVTK(target_folder_cut * "uniform_cut_$(cut_level)_data.vtu", xgrid_uni; kwargs...)
           

            ## write material map into txt files
            component_names = ["eXX","eYY","eZZ","eYZ","eXZ","eXY"]
            for c = 1 : 6
                @info "Writing strain distribution file for $(component_names[c])..."
                filename_eAB = target_folder_cut * "uniform_cut_$(cut_level)_" * component_names[c] * ".txt"
                io = open(filename_eAB, "w")
                @printf(io, "%s\n", component_names[c])
                for n = 1 : nnodes_uni
                    @printf(io, "%.6f\n",CutSolution_ϵu.entries[(c-1)*nnodes_uni+n])
                end
                close(io)
            end
            @info "Writing strain distribution file for eYZ-eYX..."
            filename_eYZ_eYX = target_folder_cut * "uniform_cut_$(cut_level)_eYZ-eYX.txt"
            io = open(filename_eYZ_eYX, "w")
            @printf(io, "eYZ-eYX\n")
            for n = 1 : nnodes_uni
                @printf(io, "%.6f\n",CutSolution_ϵu.entries[3*nnodes_uni+n] - CutSolution_ϵu.entries[5*nnodes_uni+n])
            end
            close(io)

            ## replacing NaN with 1e30 so that min/max calculation works
            replace!(CutSolution_u.entries, NaN=>1e30)
            replace!(CutSolution_∇u.entries, NaN=>1e30)
            replace!(CutSolution_ϵu.entries, NaN=>1e30)
            replace!(CutSolution_P.entries, NaN=>1e30)

            ## plot displacement, strain and polarisation on uniform cut grid
            if do_uniformcut_plots
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
                        if abs(CutSolution_ϵu.entries[j]) < 1e10
                            for k = 1 : 6
                                ϵmax[k] = max(ϵmax[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                                ϵmin[k] = min(ϵmin[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                            end
                        end
                    end
                end
                scalarplot(xgrid_uni, view(CutSolution_u.entries,1:nnodes_uni), Plotter = Plotter; flimits = (uxmin,uxmax), title = "ux on cut", fignumber = 1)
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(target_folder_cut * "uniform_cut_$(cut_level)_ux.png")
                end
                scalarplot(xgrid_uni, view(CutSolution_u.entries,nnodes_uni+1:2*nnodes_uni), Plotter = Plotter; flimits = (uymin,uymax), title = "uy on cut", fignumber = 1)
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(target_folder_cut * "uniform_cut_$(cut_level)_uy.png")
                end
                scalarplot(xgrid_uni, view(CutSolution_u.entries,2*nnodes_uni+1:3*nnodes_uni), Plotter = Plotter; flimits = (uzmin,uzmax), title = "uz on cut", fignumber = 1)
                if isdefined(Plotter,:savefig)
                    Plotter.savefig(target_folder_cut * "uniform_cut_$(cut_level)_uz.png")
                end
                if length(Solution) > 1
                    scalarplot(xgrid_uni, CutSolution_P.entries, Plotter = Plotter; flimits = (Pmin,Pmax), title = "Polarisation on cut", fignumber = 1)
                    if isdefined(Plotter,:savefig)
                        Plotter.savefig(target_folder_cut * "uniform_cut_$(cut_level)_P.png")
                    end
                end
                for k = 1 : 6
                    scalarplot(xgrid_uni, view(CutSolution_ϵu.entries,(k-1)*nnodes_uni+1:k*nnodes_uni), Plotter = Plotter; flimits = (ϵmin[k],ϵmax[k]), title = "ϵu[$k] on cut", fignumber = 1)
                    if isdefined(Plotter,:savefig)
                        Plotter.savefig(target_folder_cut * "uniform_cut_$(cut_level)_ϵ$k.png")
                    end
                end
            end

            ## revert displacing (for next cut)
            displace_mesh!(xgrid, Solution[1]; magnify = -1)
        end
    end
end


function perform_plane_cuts(target_folder_cut, Solution, plane_points, cut_levels; strain_model = NonlinearStrain3D, cut_npoints = 100, 
    only_localsearch = true, vol_cut = 16, eps_gfind = 1e-11, Plotter = nothing)

    xgrid = Solution[1].FES.xgrid

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

        x = [[plane_points[1][1],plane_points[1][2],cut_level],[plane_points[2][1],plane_points[2][2],cut_level],[plane_points[3][1],plane_points[3][2],cut_level]]

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

    ## displace grid
    displace_mesh!(xgrid, Solution[1])

    ## cut displaced grid at plane
    for l = 1 : length(cut_levels)
        cut_level = cut_levels[l]
        @info "Cutting domain at z = $cut_level with plane equation coefficients $(plane_equation_coeffs[l])"
        @time cut_grid, xgrid_uni, xtrafo!, start_cell = get_cutgrids(xgrid, plane_equation_coeffs[l]; npoints = cut_npoints, vol_cut = vol_cut)

        # plot boundary-conforming Delaunay cut mesh (suitable for FV)
       # @info "Plotting Delaunay cut mesh..."
        # gridplot(cut_grid, Plotter = Plotter, title = "Delaunay mesh of cut", fignumber = 1)

        ## interpolate data on uniform cut_grid
        @info "Interpolating data on uniform cut mesh..."
        FES2D = FESpace{H1P1{3}}(xgrid_uni)
        FES2D_∇u = FESpace{H1P1{9}}(xgrid_uni)
        FES2D_ϵ = FESpace{H1P1{6}}(xgrid_uni)
        FES2D_P = FESpace{H1P1{1}}(xgrid_uni)
        CutSolution_u = FEVector{Float64}("u (on 2D cut at z = $(cut_level))", FES2D)
        CutSolution_∇u = FEVector{Float64}("∇u (on 2D cut at z = $(cut_level))", FES2D_∇u)
        CutSolution_ϵu = FEVector{Float64}("ϵ(u) (on 2D cut at z = $(cut_level))", FES2D_ϵ)
        CutSolution_P = FEVector{Float64}("P (on 2D cut at z = $(cut_level))", FES2D_P)
        @time interpolate!(CutSolution_u[1], Solution[1]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
        @time interpolate!(CutSolution_∇u[1], Solution[1]; operator = Gradient, xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
          
        ## calculate strain from gradient interpolation on cut
        nodevals = nodevalues(CutSolution_∇u[1], Identity)
        strain = zeros(Float64,6)
        nnodes = size(nodevals,2)
        for j = 1 : nnodes
            if nodevals[1,j] == NaN
                for k = 1 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes+j] = NaN
                end
            else
                eval_strain!(strain,view(nodevals,:,j), strain_model)
                for k = 1 : 6
                    CutSolution_ϵu.entries[(k-1)*nnodes+j] = strain[k]
                end
            end
        end
        if length(Solution) > 1
            @time interpolate!(CutSolution_P[1], Solution[2]; xtrafo = xtrafo!, start_cell = start_cell, not_in_domain_value = NaN, only_localsearch = only_localsearch, eps = eps_gfind)
        end

        ## write data into csv file
        @info "Writing data into csv file..."
        writeVTK!(target_folder_cut * "cut_$(cut_level)_data.vtu", [CutSolution_u[1],CutSolution_∇u[1],CutSolution_ϵu[1],CutSolution_P[1]]; operators = [Identity, Identity, Identity, Identity])
        writeCSV!(target_folder_cut * "cut_$(cut_level)_data.txt", [CutSolution_u[1],CutSolution_∇u[1],CutSolution_ϵu[1],CutSolution_P[1]]; operators = [Identity, Identity, Identity, Identity], seperator = "\t")

        ## replacing NaN with 1e30 so that min/max calculation works
        replace!(CutSolution_u.entries, NaN=>1e30)
        replace!(CutSolution_∇u.entries, NaN=>1e30)
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
                if abs(CutSolution_ϵu.entries[j]) < 1e10
                    for k = 1 : 6
                        ϵmax[k] = max(ϵmax[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                        ϵmin[k] = min(ϵmin[k],CutSolution_ϵu[1][(k-1)*nnodes_uni+j])
                    end
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
