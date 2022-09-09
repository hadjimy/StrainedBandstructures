abstract type StrainType end
abstract type LinearStrain <: StrainType end
abstract type LinearStrain2D <: LinearStrain end
abstract type LinearStrain3D <: LinearStrain end
abstract type NonlinearStrain <: StrainType end
abstract type NonlinearStrain2D <: NonlinearStrain end
abstract type NonlinearStrain3D <: NonlinearStrain end

# allows to include other strains like Henky strain later

@inline function eval_strain!(result, input, ::Type{LinearStrain2D}; offset = 0)
    # eps_lin(u) := 1/2 (grad(u) + grad(u)^T)
    result[offset + 1] = input[1]
    result[offset + 2] = input[4]
    result[offset + 3] = input[2] + input[3]
    return nothing
end

@inline function eval_strain!(result, input, ::Type{NonlinearStrain2D}; offset = 0)
    # first add the linear part eps_lin(u) := 1/2 (grad(u) + grad(u)^T) of the strain tensor
    eval_strain!(result, input, LinearStrain2D)
    # add nonlinear part 1/2 * (grad(u)'*grad(u)) of the strain tensor
    result[offset + 1] += 1//2 * (input[1]^2 + input[3]^2)
    result[offset + 2] += 1//2 * (input[2]^2 + input[4]^2)
    result[offset + 3] += input[1]*input[2] + input[3]*input[4]
    return nothing
end

@inline function eval_strain!(result, input, ::Type{LinearStrain3D}; offset = 0)
    # eps_lin(u) := 1/2 (grad(u) + grad(u)^T)
    result[offset + 1] = input[1]
    result[offset + 2] = input[5]
    result[offset + 3] = input[9]
    result[offset + 4] = input[6] + input[8] # = 2 eps_lin(u)[2,3] = grad(u)[2,3] + grad(u)[3,2] in Voigt notation (2* 1/2 *(grad(u)[2,3] + grad(u)[3,2]))
    result[offset + 5] = input[3] + input[7] # = 2 eps_lin(u)[1,3]
    result[offset + 6] = input[2] + input[4] # = 2 eps_lin(u)[1,2]
    return nothing
end

@inline function eval_strain!(result, input, ::Type{NonlinearStrain3D}; offset = 0)
    # first add the linear part eps_lin(u) := 1/2 (grad(u) + grad(u)^T) of the strain tensor
    eval_strain!(result, input, LinearStrain3D)
    # add nonlinear part 1/2 * (grad(u)'*grad(u)) of the strain tensor
    result[offset + 1] += 1//2 * (input[1]^2 + input[4]^2 + input[7]^2)
    result[offset + 2] += 1//2 * (input[2]^2 + input[5]^2 + input[8]^2)
    result[offset + 3] += 1//2 * (input[3]^2 + input[6]^2 + input[9]^2)
    result[offset + 4] += input[3]*input[2] + input[6]*input[5] + input[9]*input[8]
    result[offset + 5] += input[1]*input[3] + input[4]*input[6] + input[7]*input[9]
    result[offset + 6] += input[1]*input[2] + input[4]*input[5] + input[7]*input[8]
    return nothing
end
