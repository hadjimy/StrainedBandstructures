
## computes the curvature (of the circle connecting points at front and back)
## and bending angle and (if given) compares the curvature to the analytical value
function compute_statistics(xgrid, Displacement::FEVectorBlock{T,Tv,Ti,FEType,APT}, scaling) where {T,Tv,Ti,FEType,APT}
    xCoordinates = xgrid[Coordinates]
    nnodes = size(xCoordinates,2)
    ncomponents = get_ncomponents(FEType)

    # find vertex number closest to (0,d/2)
    # find vertex number farthest away from (L,d/2)
    front_level = ncomponents == 2 ? [0, scaling[1]/2] : [scaling[1]/2, scaling[2]/2, 0]
    back_level = ncomponents == 2 ? [scaling[2], scaling[1]/2] : [scaling[1]/2, scaling[2]/2, scaling[3]]
    origin_point::Int = 0
    farthest_point::Int = 0
    closest_distance::Float64 = 1e30
    farthest_distance::Float64 = 1e30
    dist::Float64 = 0
    for j = 1 : nnodes
        dist = 0
        for k = 1 : ncomponents
            dist += (front_level[k] - xCoordinates[k,j])^2
        end
        if dist < closest_distance
            closest_distance = dist
            origin_point =  j
        end
        dist = 0
        for k = 1 : ncomponents
            dist += (back_level[k] - xCoordinates[k,j])^2
        end
        if dist < farthest_distance
            farthest_distance = dist
            farthest_point = j
        end
    end

    ## compute bending arc at midoint
    nodevals = nodevalues(Displacement; continuous = true)
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