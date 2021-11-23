

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

    for itet=1:size(xCellNodes,2)
        for i=1:4
            node_indices[i]=xCellNodes[i,itet]
        end
        
        @views map!(inode->plane_equation(plane_equation_coeffs,xCoordinates[:,inode]),planeq,node_indices)
        nxs=tet_x_plane!(ixcoord,ixvalues,iedges,xCoordinates,node_indices,planeq,func; tol=tol)
        # check if cutted edges are boundary edges
        for iedge = 1 : nxs
            ixbnd[iedge] = edge_on_boundary[xCellEdges[iedges[iedge],itet]]
        end
        # save triangles of cut
        pushtris(nxs,ixcoord,ixvalues,ixbnd)
    end
    all_ixcoord, all_ixfaces, all_ixvalues

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
    z = xCoordinates[end,3]
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
    
    return xgrid, R, z
end



## computes two grids: one boundary conforming Delaunay grid for finite volume methods
## and one uniform non-boundary-conforming one that comprises the cut;
## also a function is returned that transforms 3D coordinates to 2D coordinates on the cut
function get_cutgrids(xgrid, plane_equation_coeffs; npoints = 100, vol_cut = 1.0)
    @time cut_grid, R, z = plane_cut(xgrid, plane_equation_coeffs; remesh = true, vol = vol_cut)

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

    return cut_grid, xgrid_uni, xtrafo!
end