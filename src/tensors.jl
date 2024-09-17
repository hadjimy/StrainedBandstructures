"""
$(TYPEDEF)
"""
abstract type AbstractTensor{T} end

"""
    $(TYPEDEF)

Isotropic elasticity tensor with two Lame parameters λ and μ.
"""
struct IsotropicElasticityTensor{dim,T} <: AbstractTensor{T} 
    λ::T
    μ::T
    IsotropicElasticityTensor(λ,μ,dim) = new{dim, Float64}(λ,μ)
end

"""
````
function IsotropicElasticityTensor(MD::MaterialParameters)
````

Construct an isotropic elasticity tensor from the given MaterialParameters.

"""
function IsotropicElasticityTensor(MD::MaterialParameters{T, MT, MST}) where {T,MT,MST}
    # compute Lamé constants in N/(m)^2
    if haskey(d, "λ") && haskey(d, "μ")
        λ = MD("λ")
        μ = MD("μ")
    elseif haskey(d, "C11") && haskey(d, "C44")
        λ = MD("C11") - 2*MD("C44")
        μ = MD("C44")
    else
        @error "material data does not contain enough values to init isotropic elasticity tensor (needs λ, μ or C11, C44)"
    end
    return IsotropicElasticityTensor{3,T}(λ, μ)
end


"""
    $(TYPEDEF)

Custom tensor with sparse matrix representation
"""
struct SparseMatrixTensor{T} <: AbstractTensor{T}
    C::SparseMatrixCSC{T,Int}
end

function SparseMatrixTensor(M::AbstractMatrix{T}) where {T}
    return SparseMatrixTensor{T}(SparseMatrixCSC{T,Int}(M))
end


"""
    $(TYPEDEF)

Custom tensor with dense matrix representation
"""
struct DenseMatrixTensor{T} <: AbstractTensor{T}
    C::Matrix{T}
end

# TODO: Wurtzite/ZincBlende etc. could be cast into a constructor for CustomMatrixTensor
#       or we could even make subtypes with apply functions for them to skip the zeros in their evaluation


@inline function apply_tensor!(result, input, ST::IsotropicElasticityTensor{2,T}, offset = 0, input_offset = 0) where {T}
    result[offset + 1] = ST.λ * (input[input_offset + 1] + input[input_offset + 2]) + 2 * ST.μ * input[input_offset + 1]
    result[offset + 2] = ST.λ * (input[input_offset + 1] + input[input_offset + 2]) + 2 * ST.μ * input[input_offset + 2]
    result[offset + 3] = 2 * ST.μ * input[input_offset + 3]
    return nothing 
end

@inline function apply_tensor!(result, input, ST::IsotropicElasticityTensor{3,T}, offset = 0, input_offset = 0) where {T}
    result[offset + 1] = ST.λ * (input[input_offset + 1] + input[input_offset + 2] + input[input_offset + 3]) + 2 * ST.μ * input[input_offset + 1]
    result[offset + 2] = ST.λ * (input[input_offset + 1] + input[input_offset + 2] + input[input_offset + 3]) + 2 * ST.μ * input[input_offset + 2]
    result[offset + 3] = ST.λ * (input[input_offset + 1] + input[input_offset + 2] + input[input_offset + 3]) + 2 * ST.μ * input[input_offset + 3]
    result[offset + 4] = 2 * ST.μ * input[input_offset + 4]
    result[offset + 5] = 2 * ST.μ * input[input_offset + 5]
    result[offset + 6] = 2 * ST.μ * input[input_offset + 6]
    return nothing 
end

@inline function apply_tensor!(result, input, ST::SparseMatrixTensor{T}, offset = 0, input_offset = 0) where {T}
    C = ST.C
    mul!(view(result,offset+1:offset+size(C,1)), C, view(input, input_offset+1:input_offset+size(C,2)))
    return nothing 
end

@inline function apply_tensor!(result, input, ST::DenseMatrixTensor{T}, offset = 0, input_offset = 0) where {T}
    C = ST.C
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
