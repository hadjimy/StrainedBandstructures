
## computes the curvature (of the circle connecting points at front and back)
## and bending angle and (if given) compares the curvature to the analytical value
#function compute_statistics(xgrid, Displacement::FEVectorBlock{T,Tv,Ti,FEType,APT}, scaling) where {T,Tv,Ti,FEType,APT}
function compute_statistics(xgrid, Displacement, bending_axis_end_points, FEType)
    xCoordinates = xgrid[Coordinates]
    nnodes = size(xCoordinates,2)
    ncomponents = get_ncomponents(FEType)

    repair_grid!(Displacement.FES.xgrid)

    base_level = bending_axis_end_points[1] # interface point at base cross section
    top_level = bending_axis_end_points[2]  # interface point at top cross section
    origin_point::Int = 0
    farthest_point::Int = 0
    closest_distance::Float64 = 1e30
    farthest_distance::Float64 = 1e30
    dist::Float64 = 0
    for j = 1 : nnodes
        dist = 0
        for k = 1 : ncomponents
            dist += (base_level[k] - xCoordinates[k,j])^2
        end
        if dist < closest_distance
            closest_distance = dist
            origin_point =  j
        end
        dist = 0
        for k = 1 : ncomponents
            dist += (top_level[k] - xCoordinates[k,j])^2
        end
        if dist < farthest_distance
            farthest_distance = dist
            farthest_point = j
        end
    end

    ## compute bending arc at midoint
    nodevals = nodevalues(Displacement, Identity; continuous = true)
    dist_unbend = sqrt(sum((xCoordinates[:,origin_point] - xCoordinates[:,farthest_point]).^2))
    dist_bend = sqrt(sum((xCoordinates[:,origin_point] - xCoordinates[:,farthest_point] - nodevals[:,farthest_point]).^2))
    dist_farthest = sqrt(sum((nodevals[:,farthest_point]).^2))
   # if dist_bend < 10.
   #     angle_rad = (180 - dist_bend/11)/180*π
   # else
        angle_rad = acos((dist_unbend.^2+dist_bend.^2-dist_farthest.^2)/(2*dist_unbend*dist_bend))
   # end
    angle = angle_rad * 180/pi
    #R = dist_bend*sin((π-angle_rad)/2)/sin(angle_rad)
    R = dist_bend/(2*sin(angle_rad))
    curvature = 1/R
    return angle, curvature, dist_bend, xCoordinates[:,farthest_point] + nodevals[:,farthest_point]
end