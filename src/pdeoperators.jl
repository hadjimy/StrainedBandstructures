
## constructor for NonlinearForm for displacement equation
## ((1 + ∇u) ℂ(ϵ(u) - ϵ0) / (1 + α), ∇v)
##
## note: since (1 + ∇u) is not symmetric, it is tested with full gradient ∇v

struct PDEDisplacementOperator{T, ST <: StrainType, STTType, FTypeVoigtCompress <: Function, FTypeAddGradUxStress <: Function}
    STT::STTType     # stress tensor ℂ
    misfit_strain::Vector{Vector{T}}
    α::Vector{T}
    emb::Vector{T}
    uncompress_voigt!::FTypeVoigtCompress
    add_gradu_x_stress!::FTypeAddGradUxStress
    dim::Int
    cache_offset::Int
end


(op::PDEDisplacementOperator{T,ST})(result, input, item) where {T,ST} = begin
    ## input = [∇u] written as a vector of length dim^2
    ## result = (1 + ∇u) ℂϵ(u) / (1 + α) written as a vector of length dim^2 (to be multiplied by ∇v)

    ## evaluate strain into result (Voigt notation)
    eval_strain!(result, input, ST)

    ## subtract misfit strain
    @views result[1:op.dim] .-= op.misfit_strain[item[3]]

    ## multiply strain in result with isotropic stress tensor
    ## and store in input cache (all in Voigt notation)
    apply_elasticity_tensor!(input, result, op.STT[item[3]]; offset = op.cache_offset)

    ## uncompress stress in Voigt notation (input) into full matrix (result)
    op.uncompress_voigt!(result, input; offset = op.cache_offset)

    if op.emb[1] > 0.
        # add emb[1]*grad(u)*sigma(u) (in case of complicated model)
        op.add_gradu_x_stress!(result, input; factor = op.emb[1], offset = op.cache_offset)
    end

    ## multiply by -1/(1 + α) * I
    result .*= -1 / (1 .+ op.α[item[3]])

    return nothing
end

function PDEDisplacementOperator(STT, ST, misfit_strain, α, emb, dim)
    uncompress_voigt! = (dim == 2) ? uncompress_voigt2D! : uncompress_voigt3D!
    add_gradu_x_stress! = (dim == 2) ? add_gradu_x_stress2D! : add_gradu_x_stress3D!
    return PDEDisplacementOperator{Float64, ST, typeof(STT), typeof(uncompress_voigt!), typeof(add_gradu_x_stress!)}(STT,misfit_strain,α,emb,uncompress_voigt!,add_gradu_x_stress!,dim,dim^2)
end



struct PDEPolarisationOperator{T, ST <: StrainType, PETType}
    PET::PETType     # piezo-electricity tensor ℂ
    misfit_strain::Vector{Vector{T}}
    κ::Vector{T}
    F::Matrix{T}     # storage for F, okay in this way for linear operators
    dim::Int
    cache_offset::Int
end


(op::PDEPolarisationOperator{T,ST})(result, input, item) where {T,ST} = begin
    ## input = [∇u, ∇P] written as a vector of length dim^2 + dim
        # result = E eps(u) - κ * det(F).inv(F).inv(F)^(T) * grad(V_P), where F = 1 + ∇u

        polarisation_offset = op.dim^2

        ## evaluate strain (Voigt notation)
        eval_strain!(input, input, ST; offset = op.cache_offset)

        ## subtract misfit strain
        input[op.cache_offset+1:op.cache_offset+op.dim] .-= op.misfit_strain[item[3]]

        # apply piezo-electricity tensor (result = E eps(u))
        apply_piezoelectricity_tensor!(result, input, op.PET[item[3]]; input_offset = op.cache_offset)

        if (false)        
            # construct F = I + grad(u)
            F = zeros(typeof(input[1]), op.dim, op.dim)
            for i = 1 : op.dim, j = 1 : op.dim
                F[i,j] = input[op.dim*(i-1)+j]
            end
            for i = 1 : op.dim
                F[i,i] += 1
            end
            # subtract κ det(F) inv(FF^T) * ∇P
            result .-= op.κ[item[3]] * det(F) * inv(transpose(F)*F) * view(input,polarisation_offset+1:polarisation_offset+op.dim)
        else
            # subtract κ ∇P
            result .-= op.κ[item[3]] * view(input,polarisation_offset+1:polarisation_offset+op.dim)
        end

    return nothing
end

function PDEPolarisationOperator(PET, ST, misfit_strain, κ, dim)
    return PDEPolarisationOperator{Float64, ST, typeof(PET)}(PET,misfit_strain,κ,zeros(Float64, dim, dim),dim,dim+dim^2)
end




function get_displacement_operator_new(
    STT::Vector{<:ElasticityTensorType},     # stress tensor ℂ
    ST::Type{<:StrainType},     # strain type ϵ (e.g. Nonlinear3D)
    misfit_strain,              # misfit strain ϵ0 (constant scalar or vector)
    α;                          # average values (constant scalar or vector)
    displacement_id = 1,        # unknown id for displacement
    dim = 2,                    # dimension
    emb = [1],                  # embedding coefficients for complicated nonlinear features
                                # (overwrite with your array to use with embedding solver)
    regions = [0],              # regions where the operator integrates
    bonus_quadorder = 2)        # quadrature order
    
    return NonlinearForm(Gradient,                                      # operator for test function
                         [Gradient],                                    # operators for ansatz function
                         [displacement_id],                             # unknown_ids for ansatz function
                         PDEDisplacementOperator(STT, ST, misfit_strain, α, emb, dim),
                         [dim^2, dim^2, Int(dim^2+dim*(dim+1)/2)];      # argument sizes for kernel function result and input and cache
                         name = "((I + emb*∇u)C(ϵ(u)-ϵ0) : ∇v)",         # name for print-outs
                         regions = regions,                             # regions where nonlinearform intergrates
                         dependencies = "I",                            # make it item-dependent (to get the region)
                         bonus_quadorder = bonus_quadorder,             # quadrature order
                         store = false,
                         newton = true)                                 # activate Newton derivatives (false won't work here)
end

function get_displacement_operator(
    STT::ElasticityTensorType,     # stress tensor ℂ
    ST::Type{<:StrainType},     # strain type ϵ (e.g. Nonlinear3D)
    misfit_strain,              # misfit strain ϵ0 (constant scalar or vector)
    α;                          # average values (constant scalar or vector)
    displacement_id = 1,        # unknown id for displacement
    dim = 2,                    # dimension
    emb = [1],                  # embedding coefficients for complicated nonlinear features
                                # (overwrite with your array to use with embedding solver)
    regions = [0],              # regions where the operator integrates
    store = Threads.nthreads() > 1,  # separate storage for operator (allows parallel assembly, but only reasonable if nthreads() > 1)
    bonus_quadorder = 2)        # quadrature order

    cache_offset::Int = dim^2
    uncompress_voigt! = (dim == 2) ? uncompress_voigt2D! : uncompress_voigt3D!
    add_gradu_x_stress! = (dim == 2) ? add_gradu_x_stress2D! : add_gradu_x_stress3D!

    function nonlinear_operator_kernel!(result, input)
        ## input = [∇u] written as a vector of length dim^2
        ## result = (1 + ∇u) ℂ ϵ(u) / (1 + α) written as a vector of length dim^2 (to be multiplied by ∇v)

        ## evaluate strain into result (Voigt notation)
        eval_strain!(result, input, ST)

        ## subtract misfit strain
        @views result[1:dim] .-= misfit_strain

        ## multiply strain in result with isotropic stress tensor
        ## and store in input cache (all in Voigt notation)
        apply_elasticity_tensor!(input, result, STT; offset = cache_offset)

        ## uncompress stress in Voigt notation (input) into full matrix (result)
        uncompress_voigt!(result, input; offset = cache_offset)

        if emb[1] > 0.
            # add emb[1]*grad(u)*sigma(u) (in case of complicated model)
            add_gradu_x_stress!(result, input; factor = emb[1], offset = cache_offset)
        end

        ## multiply by -1/(1 + α) * I
        result .*= -1 / (1 .+ α)

        return nothing
    end
    return NonlinearForm(Gradient,                                      # operator for test function
                         [Gradient],                                    # operators for ansatz function
                         [displacement_id],                             # unknown_ids for ansatz function
                         nonlinear_operator_kernel!,                    # kernel function (above)
                         [dim^2, dim^2, Int(dim^2+dim*(dim+1)/2)];      # argument sizes for kernel function result and input and cache
                         name = "(I + emb*∇u)C(ϵ(u)-ϵ0) : ∇v) $(store ? "[stored]" : "")",         # name for print-outs
                         regions = regions,                             # regions where nonlinearform intergrates
                         bonus_quadorder = bonus_quadorder,             # quadrature order
                         store = store,
                         newton = true)                                 # activate Newton derivatives (false won't work here)
end

## assembles bilinearform (for the displacement off-diagonal block of the polarisation equation)
##      (𝔼(ϵ(u) - ϵ0), ∇W_P)
## where F = I + grad(u)
function get_polarisation_from_strain_operator(
    PET::PiezoElectricityTensorType,    # piezo electricity tensor 𝔼
    ST::Type{<:StrainType},             # strain type ϵ (e.g. Nonlinear3D)
    misfit_strain;                      # misfit strain ϵ0 (constant scalar or vector)
    dim = 2,                            # dimension
    regions = [0],                      # regions where the operator integrates
    quadorder = 0)                      # (extra) quadrature order

    size_strain = Int(dim*(dim+1)/2)
    strain = zeros(Float64,size_strain)

    function closure(result, input)
        # input = grad(u) written as a vector
        # result = E eps(u)

        ## evaluate strain into result (Voigt notation)
        eval_strain!(strain, input, ST)

        ## subtract misfit strain
        strain .-= misfit_strain

        # apply piezo-electricity tensor
        apply_piezoelectricity_tensor!(result, strain, PET)

        return nothing
    end
    polarisation_action = Action(closure, [dim,dim^2]; dependencies = "", bonus_quadorder = quadorder)
    return BilinearForm( [Gradient,Gradient],
                         polarisation_action;
                         name = "-E (ϵ(u) - ϵ0) ⋅ ∇(Q)",
                         factor = -1,
                         regions = regions,
                         apply_action_to = [2])
end

## assembles trilinearform (for the diagonal block of the polarisation equation)
##      (det(F) inv(FF^(T)) * ∇V_P, ∇W_P)
## where F = I + ∇u (so assembly needs ∇u and ∇V_P)
function get_polarisation_laplacian_operator(;
    displacement_id = 1,        # unknown id for displacement
    dim = 2,                    # dimension
    κ = 1,                      # dielectric constant x relative permittivity
    regions = [0],              # regions where the operator integrates
    quadorder = 0)              # (extra) quadrature order

    F = zeros(Float64, dim, dim) # storage for F, okay in this way for linear operators
    function closure(result, input)
        # input = [∇u,∇V_P] written as a vector
        #       = [ 1:9   ,   10:12   ]

        # construct F = I + grad(u)
        if (true)
            for i = 1 : dim, j = 1 : dim
                F[i,j] = input[3*(i-1)+j]
            end
            for i = 1 : dim
                F[i,i] += 1
            end
            # compute result (allocations !!! idea: exploit FF^T = 2*ϵ(u) ?)
            result .= det(F)*inv(transpose(F)*F) * view(input,10:12)
        else
            result .= view(input,10:12)
        end

        return nothing
    end
    potential_action = Action(closure, [dim, dim^2+dim]; dependencies = "", bonus_quadorder = quadorder)
    return BilinearForm( [Gradient,Gradient],
                         [Gradient],
                         [displacement_id],
                          potential_action;
                          name = "-κ det(F) inv(FF^T) ∇(V_P) ⋅ ∇(Q)",
                          factor = κ,
                          regions = regions)
end

#= constructs ItemIntegrator that evaluates the stored energy
  
   E(u) = 1/8 ∫ det(M) (F^TF - I) : ℂ (F^TF - I)
  
   where F := (1 + ∇u) inv(M)
   and   M := 1 + α I

   # todo: add polarisation term, seems to work in 2D, but not in 3D
   # wish: would be nice to generate PDE operators by autodiff from given energy
=#
function get_energy_integrator(
    STT::Array{<:ElasticityTensorType,1},      # stress tensor ℂ
    ST::Type{<:StrainType},     # strain type (e.g. Nonlinear3D)
    α;                          # average values (constant scalar or vector)
    dim = 2,                    # dimension
    regions = 1:3,                    # dimension
    quadorder = 4)              # quadrature order

    ## prepare energy calculation
    size_epsu = Int(dim*(dim+1)/2)
    temp = zeros(Float64,size_epsu)
    temp2 = zeros(Float64,size_epsu)

    function energy_kernel(result, input, region)
        # input = [∇u] as vector

        # compute Voigt vector of (F^TF - I) 
        # = inv(M)^2*(∇u + ∇u' + ∇u'∇u + 1) - 1
        # = inv(M)^2*(2*ϵ(u) + 1) - 1
        # inv(M) = 1/(1+a)
        eval_strain!(temp, input, ST)
        temp .*= 2
        @views temp[1:dim] .+= 1
        temp ./= (1 + α[region])^2
        @views temp[1:dim] .-= 1

        ## multiply with isotropic stress tensor
        fill!(temp2,0)
        apply_elasticity_tensor!(temp2, temp, STT[region])

        ## compute energy
        result[1] = 0
        for j = 1 : size_epsu
            result[1] += temp2[j] * temp[j]
        end

        ## multiply det(M)/8
        result[1] *= (1 + α[region])^2 / 8

        return nothing  
    end
    energy_action = Action{Float64}(energy_kernel, [1,dim^2]; name = "stored energy", dependencies = "R", bonus_quadorder = quadorder)
    return ItemIntegrator([Gradient], energy_action; regions = regions)
end



## WIP
## constructor for NonlinearForm for polarisation equation
## (very slow, better use the linear ones above)
function get_polarisation_operator(
    PET::PiezoElectricityTensorType,      # piezo electricity tensor 𝔼
    ST::Type{<:StrainType},     # strain type ϵ (e.g. Nonlinear3D)
    misfit_strain;              # misfit strain ϵ0 (constant scalar or vector)
    displacement_id = 1,        # unknown id for displacement u
    polarisation_id = 2,        # unknown id for polarisation potential V_P
    dim = 2,                    # dimension
    κ = 1,                      # dielectric constant (times relative permittivity)
    emb = [1],                  # embedding coefficients for complicated nonlinear features
                                # (overwrite with your array to use with embedding solver)
    regions = [0],              # regions where the operator integrates
    quadorder = 2)              # quadrature order

    polarisation_offset = dim^2
    cache_offset::Int = dim^2 + dim
    F = zeros(Real,dim^2,dim^2)

    function nonlinear_operator_kernel!(result, input)
        ## input = [grad(u), grad(V_P)] written as a vector
        # result = E eps(u) - κ * det(F).inv(F).inv(F)^(T) * grad(V_P), where F = 1 + ∇u

        ## evaluate strain into result (Voigt notation)
        eval_strain!(input, input, ST; offset = cache_offset)

        ## subtract misfit strain
        input[cache_offset+1:cache_offset+dim] .-= misfit_strain

        # apply piezo-electricity tensor
        apply_piezoelectricity_tensor!(result, input, PET; input_offset = cache_offset)


        # TODO : repair this such that it works for DualNumbers
        ## build F = 1 + ∇u
        #for i = 1 : dim, j = 1 : dim
        #    F[i,j] = input[dim*(i-1)+j]
        #end
        #for i = 1 : dim
        #    F[i,i] += 1
        #end

        # subtract κ det(F) inv(FF^T) * grad(V_P)
        #result .-= κ * det(F) * inv(transpose(F)*F) * view(input,polarisation_offset+1:polarisation_offset+dim)
        result .-= κ * view(input,polarisation_offset+1:polarisation_offset+dim)

        return nothing
    end
    return NonlinearForm(Gradient,                                  # operator for test function
                         [Gradient, Gradient],                      # operators for ansatz function
                         [displacement_id, polarisation_id],        # unknown_ids for ansatz function
                         nonlinear_operator_kernel!,                # kernel function (above)
                         [dim, dim^2+dim, Int(dim^2+dim+dim*(dim+1)/2)];          # argument sizes for kernel function result and input and cache
                         name = "(E ϵ(u) - κ det(F) inv(FF^T) ∇V_P) : ∇W_P",     # name for print-outs
                         regions = regions,                         # regions where nonlinearform intergrates
                         bonus_quadorder = quadorder,                     # quadrature order
                         newton = true)                             # activate Newton derivatives (false won't work here)
end

# Voigt notation compresses entries of symmetric matrices
#
# 2D Voigt(S) [1,2,3] --> [1 3
#                          3 2]
#
@inline function uncompress_voigt2D!(result, input; offset = 0)
    result[1] = input[offset+1]
    result[2] = input[offset+3]
    result[3] = input[offset+3]
    result[4] = input[offset+2]
    return nothing
end

# Voigt notation compresses entries of symmetric matrices
#
#                               [1 6 5
# 3D Voigt(S) [1,2,3,4,5,6] -->  6 2 4
#                                5 4 3]
@inline function uncompress_voigt3D!(result, input; offset = 0)
    result[1] = input[offset+1]
    result[2] = input[offset+6]
    result[3] = input[offset+5]
    result[4] = input[offset+6]
    result[5] = input[offset+2]
    result[6] = input[offset+4]
    result[7] = input[offset+5]
    result[8] = input[offset+4]
    result[9] = input[offset+3]
    return nothing
end

# some helper functions to compute 
#    result .+= ∇u σ(u) (matrix-matrix product)
# input = [∇u, σ(u)]
@inline function add_gradu_x_stress2D!(result, input; offset = 0, factor = 1)
    result[1] += factor*(input[offset+1]*input[1] + input[offset+3]*input[2])
    result[2] += factor*(input[offset+3]*input[1] + input[offset+2]*input[2])
    result[3] += factor*(input[offset+1]*input[3] + input[offset+3]*input[4])
    result[4] += factor*(input[offset+3]*input[3] + input[offset+2]*input[4])
    return nothing
end
@inline function add_gradu_x_stress3D!(result, input; offset = 0, factor = 1)
    result[1] += factor * (input[offset+1]*input[1] + input[offset+6]*input[2] + input[offset+5]*input[3])
    result[2] += factor * (input[offset+6]*input[1] + input[offset+2]*input[2] + input[offset+4]*input[3])
    result[3] += factor * (input[offset+5]*input[1] + input[offset+4]*input[2] + input[offset+3]*input[3])
    result[4] += factor * (input[offset+1]*input[4] + input[offset+6]*input[5] * input[offset+5]*input[6])
    result[5] += factor * (input[offset+6]*input[4] + input[offset+2]*input[5] * input[offset+4]*input[6])
    result[6] += factor * (input[offset+5]*input[4] + input[offset+4]*input[5] * input[offset+3]*input[6])
    result[7] += factor * (input[offset+1]*input[7] + input[offset+6]*input[8] * input[offset+5]*input[9])
    result[8] += factor * (input[offset+6]*input[7] + input[offset+2]*input[8] * input[offset+4]*input[9])
    result[9] += factor * (input[offset+5]*input[7] + input[offset+4]*input[8] * input[offset+3]*input[9])
    return nothing
end
