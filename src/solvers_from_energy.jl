mutable struct EnergyOperator{MT, STTType}
    M::Array{Matrix{MT},1}
    detM::Array{MT,1}
    invM::Array{Matrix{MT},1}
    STT::STTType
    eval_W!::Function
	eval_∂FW!::Function
end

function update_M!(EO::EnergyOperator, M)
    EO.M .= M
    EO.detM .= [det(m) for m in M]
    EO.invM .= [inv(m) for m in M]
    EO.eval_W!, EO.eval_∂FW! = ∂FW!(EO.M, EO.invM, EO.detM, EO.STT)
end

function EnergyOperator(M, STT)
    detM = [det(m) for m in M]
    invM = [inv(m) for m in M]
    eval_W!, eval_∂FW! = ∂FW!(M, invM, detM, STT)
    return EnergyOperator(M, detM, invM, STT, eval_W!, eval_∂FW!)
end


## derivative of energy (by ForwardDiff)
function ∂FW!(M, invM, detM, STT)
    Dresult = nothing
	cfg = nothing
	result_dummy = nothing
    Fel = nothing
    ϵ = nothing
    Cϵ = nothing

    function Wgeneral!(result, _ϵ, _Cϵ, _F, region)
        # compute strain ϵ := (F'F - I)/2
        compute_strain!(_ϵ, _F)

        # apply elasticity tensor Cϵ = C : ϵ
        apply_elasticity_tensor!(_Cϵ, _ϵ, STT[region])

        # compute energy
        result[1] = dot(_Cϵ, _ϵ) / 2
    end

    function _W!(result, Fel, ϵ, Cϵ, F, qpinfo)
        ## extract region number
        region = qpinfo.region

        ## compute elastic strain F_el = F*inv(M)
        multiply_matrices_as_vectors!(Fel, F, invM[region])

        ## evaluate general energy
        Wgeneral!(result, ϵ, Cϵ, Fel, region)

        ## multiply det(M)
        result .*= detM[region]
    end

    function W!(result, _F, qpinfo)
        if Fel == nothing
            Fel = zeros(eltype(_F), 9)
            ϵ = zeros(eltype(_F), 6)
            Cϵ = zeros(eltype(_F), 6)
        end
        _W!(result, Fel, ϵ, Cϵ, _F, qpinfo)
    end

	energy(qpinfo) = (a, b) -> W!(a, b, qpinfo)
	
    function _closure(result, ∇u, qpinfo, result_dummy, Dresult, cfg)
		Dresult = ForwardDiff.vector_mode_jacobian!(Dresult, energy(qpinfo), result_dummy, ∇u, cfg)
		result .= view(DiffResults.jacobian(Dresult), :)
    end

	function closure(result, ∇u, qpinfo)

		if Dresult === nothing
			## first initialization of DResult when type of F is known
			result_dummy = zeros(eltype(∇u), 1)
			Dresult = DiffResults.JacobianResult(result_dummy, ∇u)
			cfg = ForwardDiff.JacobianConfig(energy(qpinfo), result_dummy, ∇u, ForwardDiff.Chunk{length(∇u)}())
		end

        ## compute F := I + ∇u from ∇u (below ∇u acts as F)
        ∇u[1] += 1 
        ∇u[5] += 1
        ∇u[9] += 1

        ## we can take the derivative with respect to ∇u, since F = I + ∇u and therefore ∂F = ∂F(∇u)
        _closure(result, ∇u, qpinfo, result_dummy, Dresult, cfg)
		return nothing
	end

    return W!, closure
end


function multiply_matrices_as_vectors!(A,B,C) # A = B*C
    A[1] = B[1] * C[1] + B[2] * C[4] + B[3] * C[7]
    A[2] = B[1] * C[2] + B[2] * C[5] + B[3] * C[8]
    A[3] = B[1] * C[3] + B[2] * C[6] + B[3] * C[9]
    A[4] = B[4] * C[1] + B[5] * C[4] + B[6] * C[7]
    A[5] = B[4] * C[2] + B[5] * C[5] + B[6] * C[8]
    A[6] = B[4] * C[3] + B[5] * C[6] + B[6] * C[9]
    A[7] = B[7] * C[1] + B[8] * C[4] + B[9] * C[7]
    A[8] = B[7] * C[2] + B[8] * C[5] + B[9] * C[8]
    A[9] = B[7] * C[3] + B[8] * C[6] + B[9] * C[9]
end

function compute_strain!(ϵ, F) # ϵ = (F'*F - I)/2
    ## as a Voigt vector
    ϵ[1] = F[1] * F[1] + F[4] * F[4] + F[7] * F[7] - 1
    ϵ[2] = F[2] * F[2] + F[5] * F[5] + F[8] * F[8] - 1
    ϵ[3] = F[3] * F[3] + F[6] * F[6] + F[9] * F[9] - 1
    ϵ[4] = 2*(F[2] * F[3] + F[5] * F[6] + F[8] * F[9])
    ϵ[5] = 2*(F[1] * F[3] + F[4] * F[6] + F[7] * F[9])
    ϵ[6] = 2*(F[2] * F[1] + F[5] * F[4] + F[8] * F[7])
    ϵ ./= 2
end
