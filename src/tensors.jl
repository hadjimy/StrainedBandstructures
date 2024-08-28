"""
$(TYPEDEF)
"""
abstract type ElasticityTensorType{T} end

"""
    $(TYPEDEF)

Isotropic elasticity tensor with two Lame parameters λ and μ.
"""
struct IsotropicElasticityTensor{dim,T} <: ElasticityTensorType{T} 
    λ::T
    μ::T
    IsotropicElasticityTensor(λ,μ,dim) = new{dim, Float64}(λ,μ)
end


"""
    $(TYPEDEF)

Custom elasticity tensor that takes any matrix (e.g. 3x3 in 2D and 6x6 in 3D),
that encodes the mapping of the strain (in Voigt notation) to the stress tensor (in Voigt notation)
"""
struct CustomMatrixElasticityTensor{T} <: ElasticityTensorType{T}
    C::Matrix{T}
end
# TODO: Wurtzite/ZincBlende etc. could be cast into a constructor for CustomMatrixTensor
#       or we could even make subtypes with apply functions for them to skip the zeros in their evaluation


@inline function apply_elasticity_tensor!(result, input, ST::IsotropicElasticityTensor{2,T}; offset = 0) where {T}
    result[offset + 1] = ST.λ * (input[1] + input[2]) + 2 * ST.μ * input[1]
    result[offset + 2] = ST.λ * (input[1] + input[2]) + 2 * ST.μ * input[2]
    result[offset + 3] = 2 * ST.μ * input[3]
    return nothing 
end

@inline function apply_elasticity_tensor!(result, input, ST::IsotropicElasticityTensor{3,T}; offset = 0) where {T}
    result[offset + 1] = ST.λ * (input[1] + input[2] + input[3]) + 2 * ST.μ * input[1]
    result[offset + 2] = ST.λ * (input[1] + input[2] + input[3]) + 2 * ST.μ * input[2]
    result[offset + 3] = ST.λ * (input[1] + input[2] + input[3]) + 2 * ST.μ * input[3]
    result[offset + 4] = 2 * ST.μ * input[4]
    result[offset + 5] = 2 * ST.μ * input[5]
    result[offset + 6] = 2 * ST.μ * input[6]
    return nothing 
end

@inline function apply_elasticity_tensor!(result, input, ST::CustomMatrixElasticityTensor{T}; offset = 0) where {T}
    C::Matrix{T} = ST.C
    for k = 1 : size(C,1)
        result[offset + k] = 0
        for j = 1 : size(C,2)
            if C[k,j] != 0
                result[offset + k] += C[k,j] * input[j]
            end
        end
    end
    return nothing 
end



"""
    $(TYPEDEF)

"""
abstract type PiezoElectricityTensorType{T} end

"""
    $(TYPEDEF)

Custom piezo-elecitrcity tensor that takes any matrix (e.g. 3x2 in 2D and 6x3 in 3D),
that encodes the mapping of the stress (in Voigt notation) to the polarization.
"""
struct CustomMatrixPiezoElectricityTensor{T} <: PiezoElectricityTensorType{T}
    C::Matrix{T}
end
# Wurtzite/ZincBlende etc. could be cast into a constructor for CustomMatrixTensor

@inline function apply_piezoelectricity_tensor!(result, input, PET::CustomMatrixPiezoElectricityTensor{T}; offset = 0, input_offset = 0) where {T}
    C::Matrix{T} = PET.C
    for k = 1 : size(C,1)
        result[offset + k] = 0
        for j = 1 : size(C,2)
            if C[k,j] != 0
                result[offset+k] += C[k,j] * input[input_offset + j]
            end
        end
    end
    return nothing 
end