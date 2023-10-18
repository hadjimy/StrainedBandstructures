
## solves the problem by parameter embedding
## idea: multiply difficult parts of the operators by an embedding parameter
##       solver starts by setting these parameters to zero and then uses this solution as initial guesses for
##       reassembled problem with increased parameters, repeated for nsteps until desired parameters are reached
function solve_by_embedding(
            Problem,                                        # problem description
            xgrid::ExtendableGrid{Tv,Ti},                   # grid
            emb_params;                                     # embedding parameters (operators must depend on them!)
            FETypes = [H1P1{size(xgrid[Coordinates],1)}],   # FETypes (default: P1)
            linsolver = "UMFPACK",                          # change solver (e.g. "MKLPARDISO", "UMFPACK", or any ExtendableSparse.LUFactorization)
            nsteps = ones(Int,length(FETypes)),             # number of embedding steps (parameters are scaled by nsteps equidistant steps within 0:1)
            subiterations = [1:length(FETypes)],            # maximal iterations in each embedding step
            damping = 0,                                    # damping in Newton iteration
            target_residual = 1e-12*ones(Int,length(FETypes)),
            maxiterations = 20*ones(Int,length(FETypes))) where {Tv,Ti}

    @info "solver = $linsolver"
    ## discretise the problem
    ## create finite element space (FESpace) and solution vector (FEVector)
    ## generate FESpace and FEVector for discretisation
    FES = Array{FESpace{Tv,Ti},1}(undef, length(FETypes))
    for j = 1 : length(FETypes)
        FES[j] = FESpace{FETypes[j]}(xgrid)
    end
    Solution = FEVector(FES)

    ## prepare parameter embedding
    emb_params_target = deepcopy(emb_params)
    if sum(emb_params_target) == 0.
        nsteps .= 1
        @warn "all emb_params_target == 0, reducing nsteps to 1"
    end

    residual::Float64 = 0
    for s = 1 : length(subiterations)
        for j = 1 : nsteps[s]
            emb_params .= nsteps[s] == 1 ? emb_params_target : (j-1)/(nsteps[s]-1) .* emb_params_target

            ## solve by GradientRobustMultiPhysics standard fixed-point solver
            println("Solving problem with parameter emb_params = $emb_params (embedding step $j/$(nsteps[s]))...")
            residual = solve!(Solution, Problem; subiterations = subiterations[s], show_statistics = true, damping = damping, linsolver = linsolver, maxiterations = maxiterations[s], target_residual = target_residual[s])
        end
    end

    return Solution, residual
end

function solve_by_embedding!(
    Solution,
    Problem,                                        # problem description
    xgrid::ExtendableGrid{Tv,Ti},                   # grid
    emb_params;                                     # embedding parameters (operators must depend on them!)
    FETypes = [H1P1{size(xgrid[Coordinates],1)}],   # FETypes (default: P1)
    linsolver = "UMFPACK",                          # change solver (e.g. "MKLPARDISO", "UMFPACK", or any ExtendableSparse.LUFactorization)
    nsteps = ones(Int,length(FETypes)),             # number of embedding steps (parameters are scaled by nsteps equidistant steps within 0:1)
    subiterations = [1:length(FETypes)],            # maximal iterations in each embedding step
    damping = 0,                                    # damping in Newton iteration
    target_residual = 1e-12*ones(Int,length(FETypes)),
    maxiterations = 20*ones(Int,length(FETypes))) where {Tv,Ti}

@info "solver = $linsolver"

## prepare parameter embedding
emb_params_target = deepcopy(emb_params)
if sum(emb_params_target) == 0.
nsteps .= 1
@warn "all emb_params_target == 0, reducing nsteps to 1"
end

residual::Float64 = 0
for s = 1 : length(subiterations)
for j = 1 : nsteps[s]
    emb_params .= nsteps[s] == 1 ? emb_params_target : (j-1)/(nsteps[s]-1) .* emb_params_target

    ## solve by GradientRobustMultiPhysics standard fixed-point solver
    println("Solving problem with parameter emb_params = $emb_params (embedding step $j/$(nsteps[s]))...")
    residual = solve!(Solution, Problem; subiterations = subiterations[s], show_statistics = true, damping = damping, linsolver = linsolver, maxiterations = maxiterations[s], target_residual = target_residual[s])
end
end

return Solution, residual
end



## solves the problem by criterion-dependent damping
## idea: damp each Newton iteration if descent criterions are not satisfied
##       two criterions are available:
##          - use_energy_decrease_criterion (requires that energy decreases)
##          - use_residual_decrease_criterion (requires that nonlinear residual ≈ derivative of energy decreases)
## 
## (first criterion requires an energy integrator that should match the implemented operators)
function solve_by_damping(
            Problem,                                            # problem description
            xgrid::ExtendableGrid{Tv,Ti},                       # grid
            linsolver = "UMFPACK",                              # change solver (e.g. "MKLPARDISO", "UMFPACK", or any ExtendableSparse.LUFactorization)
            EIntegrator = nothing;                              # energy integrator
            FETypes = [H1P1{size(xgrid[Coordinates],1)}],       # FETypes (default: P1)
            target_residual = 1e-12,                            # target residual
            use_energy_decrease_criterion = true,               # damping step is accepted if energy decreases
            use_residual_decrease_criterion = false,            # damping step is accepted if nonlinear residual decreases
            use_or = true,                                      # it is enough when one of the criterions above is satisfied
            damping_step = 0.05,                                # checks damping values within the range 0:damping_step:1
            maxiterations = 100) where {Tv,Ti}                  # maximal iterations

    ## generate solution vector
    FES = Array{FESpace{Tv,Ti},1}(undef, length(FETypes))
    for j = 1 : length(FETypes)
        FES[j] = FESpace{FETypes[j]}(xgrid)
    end
    Solution = FEVector(FES)
    subiterations = [1:length(FETypes)]

    #############################
    ### DAMPING CONFIGURATION ###
    #############################

    TestVector = deepcopy(Solution)
    min_energy::Array{Float64,1} = [1e30,1e30]
    energy::Array{Float64,1} = [0.0,0.0]
    b = FEVector{Float64}("DI",Solution[1].FES)
    damping_guesses = [1.0]
    append!(damping_guesses,0.0:damping_step:0.99)

    @assert use_energy_decrease_criterion || EIntegrator === nothing "need an energy integrator"
    
    function get_damping(old_iterate,newton_iterate, fixed_dofs)

        ## this function tries damping factors in 0.1 steps and
        ## takes the smallest one that satsfies criterions

        min_energy .= [1e30,1e30]

        println("                    DAMPING |    ENERGY    NL-RESIDUAL | ACCEPTED?")
        for damping::Float64 in damping_guesses
            for j = 1 : length(TestVector[1])
                TestVector[1][j] = damping * old_iterate[1][j] + (1 - damping) * newton_iterate[1][j]
            end

            if use_energy_decrease_criterion
                energy[1] = evaluate(EIntegrator, TestVector[1])
            end
            if use_residual_decrease_criterion # todo : make more general
                fill!(b.entries,0)
                for o = 1 : length(Problem.RHSOperators[1])
                    eval_assemble!(b[1],Problem.RHSOperators[1][o],TestVector)
                end
                for dof in fixed_dofs
                    if dof <= length(b.entries)
                        b.entries[dof] = 0
                    end
                end
                energy[2] = sqrt(sum(b.entries.^2))
            end
            if damping == 1.0 # save values to beat
                min_energy[1] = energy[1]
                min_energy[2] = energy[2]
                @printf("                      %.2f  | %.6e  %.4e | (values to beat)\n", damping, energy[1], energy[2])
            else
                criterion1_satisfied = energy[1] < min_energy[1]
                criterion2_satisfied = energy[2] < min_energy[2]
                if !use_energy_decrease_criterion
                    criterion1_satisfied = true
                end
                if !use_residual_decrease_criterion
                    criterion2_satisfied = true
                end
                @printf("                      %.2f  | %.6e  %.4e | %s %s\n", damping, energy[1], energy[2], criterion1_satisfied, criterion2_satisfied)
                if use_or && criterion1_satisfied && use_energy_decrease_criterion
                    criterion2_satisfied = true
                end
                if use_or && criterion2_satisfied && use_residual_decrease_criterion
                    criterion1_satisfied = true
                end
                if criterion1_satisfied * criterion2_satisfied
                    return damping
                end
            end
        end

        @warn "could not find a damping factor that satisfies criterions, taking damping = 0"
        return 0
    end
    
    ## solve
    println("Solving by damping...")
    residual = solve!(Solution, Problem; linsolver = linsolver, subiterations = subiterations, target_residual = target_residual, maxiterations = maxiterations, damping = get_damping)

    return Solution, residual
end



function solve_lowlevel(
    xgrid::ExtendableGrid,                          # grid
    boundary_conditions,                            # Array of BoundaryData for each unknown
    global_constraints,                             # Array of lgobal constraints (like periodic dofs, fixed integral means)
    displacement_operator,                          # Operator for displacement
    polarisation_operator,                          # Operator for polarisation
    emb_params;                                     # embedding parameters (operators must depend on them!)
    FETypes = [H1P1{size(xgrid[Coordinates],1)}],   # FETypes (default: P1)
    linsolver = ExtendableSparse.MKLPardiso,
    nsteps = ones(Int,length(FETypes)),             # number of embedding steps (parameters are scaled by nsteps equidistant steps within 0:1)
    damping = 0,                                    # damping in Newton iteration
    target_residual = 1e-12*ones(Int,length(FETypes)),
    solve_polarisation = false,                     # also solve for polarisation ?
    coupled = false,                                # solve coupled system in each nonlinear iteration ? (not needed as long as polarisation does not influence displacement)
    bonus_quadorder = 2,
    use_sparsity_detection = false,                 # sparsity detection for jacobians
    maxiterations = 20*ones(Int,length(FETypes)),
    fixed_penalty = 1e60
    )

    ## generate FES spaces
    Tv = Float64
    Ti = Int32
    FES = Array{FESpace{Tv,Ti},1}(undef, 0)
    push!(FES, FESpace{FETypes[1]}(xgrid))

    solve_polarisation = (length(FETypes) == 2) * solve_polarisation
    coupled = coupled * solve_polarisation

    if solve_polarisation
        push!(FES, FESpace{FETypes[2]}(xgrid))
    end

    ## generate FEVectors and FEMatrix
    SolutionOld = FEVector(FES)     # storage for previous solution
    Solution::FEVector{Float64} = FEVector(FES)        # storage for current and final solution
    if coupled
        S = FEMatrix(FES)       # system matrix
        rhs = FEVector(FES)     # system rhs
    else
        S = FEMatrix(FES[1])       # system matrix
        rhs = FEVector(FES[1])     # system rhs
    end

    @info "Starting low level solver...
            ndofs = $(length(Solution.entries))
            polarisation = $solve_polarisation (coupled = $coupled)"

    ## generate quadrature rule
    EG = xgrid[UniqueCellGeometries][1]
    order = get_polynomialorder(FETypes[1], EG)
    qf = QuadratureRule{Float64, EG}(2*(order-1) + bonus_quadorder)
    weights::Vector{Float64} = qf.w
    nweights::Int = length(weights)
    cellvolumes::Array{Float64,1} = xgrid[CellVolumes]
    cellregions::Array{Int,1} = xgrid[CellRegions]

    ## prepare dofmaps and FE evaluators
    ndofs_u::Int = get_ndofs(ON_CELLS, FETypes[1], EG)      # displacement
    ∇u::FEEvaluator{Float64} = FEEvaluator(FES[1], Gradient, qf)
    vals_∇u::Array{Float64,3} = ∇u.cvals
    celldofs_u::Matrix{Int32} = Matrix(FES[1][GradientRobustMultiPhysics.CellDofs])
    indices_∇u = 1:size(vals_∇u, 1)

    if solve_polarisation # polarisation
        ndofs_P::Int = get_ndofs(ON_CELLS, FETypes[2], EG)  
        ∇P::FEEvaluator{Float64} = FEEvaluator(FES[2], Gradient, qf)
        vals_∇P::Array{Float64,3} = ∇P.cvals
        celldofs_P::Matrix{Int32} = Matrix(FES[2][GradientRobustMultiPhysics.CellDofs])
        indices_∇P = 1:size(vals_∇P, 1)
    end

    ## local matrix and vector structures
    AUU::Array{Float64,2} = zeros(Float64, ndofs_u, ndofs_u)
    bU::Array{Float64,1} = zeros(Float64, ndofs_u)
    SE::ExtendableSparseMatrix{Float64,Int64} = S.entries
    SolutionE::Array{Float64,1} = Solution.entries

    offsets::Array{Int,1} = [0, FES[1].ndofs]
    if solve_polarisation
        push!(offsets, FES[1].ndofs + FES[2].ndofs)
    end

    ## prepare automatic differentation of displacement operator
    dim = displacement_operator.dim
    if solve_polarisation
        @assert dim == polarisation_operator.dim
    end
    argsizes_displacement = [dim^2,dim^2,Int(dim^2+dim*(dim+1)/2)]
    argsizes_polarisation = [dim,dim+dim^2,Int(dim+dim^2+dim*(dim+1)/2)]
    item_information::Array{Ti,1} = ones(Ti,3)
    displacement_operator4region(r) = (result,input) -> displacement_operator(result,input,r)
    polarisation_operator4region(r) = (result,input) -> polarisation_operator(result,input,r)
    result_displacement::Array{Float64,1} = Vector{Float64}(undef,argsizes_displacement[1])
    input_displacement::Array{Float64,1} = Vector{Float64}(undef,argsizes_displacement[3])
    Dresult_displacement = DiffResults.JacobianResult(result_displacement,input_displacement)
    result_polarisation::Array{Float64,1} = Vector{Float64}(undef,argsizes_polarisation[1])
    input_polarisation::Array{Float64,1} = Vector{Float64}(undef,argsizes_polarisation[3])
    Dresult_polarisation = DiffResults.JacobianResult(result_polarisation,input_polarisation)

    if (use_sparsity_detection)
        sparsity_pattern = Symbolics.jacobian_sparsity(displacement_operator4region(item_information),result_displacement,input_displacement)
        jac = Float64.(sparse(sparsity_pattern))
        colors = matrix_colors(jac) 
        cfg_displacement = ForwardColorJacCache(displacement_operator4region(item_information),input_displacement,nothing;
                dx = nothing,
                colorvec = colors,
                sparsity = nothing)
    else
        cfg_displacement = ForwardDiff.JacobianConfig(displacement_operator4region(item_information), result_displacement, input_displacement, ForwardDiff.Chunk{argsizes_displacement[3]}())
        cfg_polarisation = ForwardDiff.JacobianConfig(polarisation_operator4region(item_information), result_polarisation, input_polarisation, ForwardDiff.Chunk{argsizes_polarisation[3]}())
    end
    jac_displacement::Array{Float64,2} = DiffResults.jacobian(Dresult_displacement)
    value_displacement::Array{Float64,1} = DiffResults.value(Dresult_displacement)
    jac_polarisation::Array{Float64,2} = DiffResults.jacobian(Dresult_polarisation)
    value_polarisation::Array{Float64,1} = DiffResults.value(Dresult_polarisation)
    
    ## ASSEMBLY LOOP
    fixed_dofs = nothing
    fixed_dofsP = nothing
    eval_∇u::Array{Float64,1} = zeros(Float64, argsizes_displacement[3])
    eval_∇P::Array{Float64,1} = zeros(Float64, argsizes_polarisation[3])
    jac_displacement_view = view(jac_displacement,indices_∇u,indices_∇u)
    ncells::Int = num_cells(xgrid)
    tempV::Array{Float64,1} = zeros(Float64, argsizes_displacement[1])
    tempP::Array{Float64,1} = zeros(Float64, argsizes_polarisation[1])

    function update_displacement(; coupled = false)

        fill!(S[1,1],0)
        fill!(rhs[1],0)
        if coupled
            fill!(S[1,2],0)
        end

        for cell = 1 : ncells

            ## update FE basis evaluators
            item_information[1] = cell
            item_information[2] = cell
            item_information[3] = cellregions[cell]
            disp_op = displacement_operator4region(item_information)
            ∇u.citem[] = cell
            update_basis!(∇u) 
            dofs_u = view(celldofs_u, : , cell)

            for qp = 1 : nweights
                ## evaluate ∇u in current solution at current quadrature point
                fill!(eval_∇u, 0)
                for j = 1 : ndofs_u, k in indices_∇u
                    eval_∇u[k] += SolutionE[dofs_u[j]] * vals_∇u[k,j,qp] 
                end

                ## evaluate jacobian of displacement tensor
                if use_sparsity_detection
                    forwarddiff_color_jacobian!(jac_displacement, disp_op, eval_∇u, cfg_displacement)
                else
                    ForwardDiff.vector_mode_jacobian!(Dresult_displacement,disp_op,result_displacement,eval_∇u,cfg_displacement)
                end

                # update matrix
                for j = 1 : ndofs_u
                    # multiply with jacobian
                    mul!(tempV, jac_displacement_view, view(vals_∇u,:,j,qp))

                    # multiply test function operator evaluation
                    for k = 1 : ndofs_u
                        AUU[j,k] += dot(tempV, view(vals_∇u,:,k,qp)) * weights[qp]
                    end
                end 

                # update rhs
                mul!(tempV, jac_displacement_view, view(eval_∇u, indices_∇u))
                disp_op(value_displacement, eval_∇u)
                tempV .-= value_displacement
                for j = 1 : ndofs_u
                    bU[j] += dot(tempV, view(vals_∇u,:,j,qp)) * weights[qp]
                end
            end

            # write into full matrix
            for j = 1 : ndofs_u
                for k = 1 : ndofs_u
                    if abs(AUU[j,k]) > 0
                        rawupdateindex!(SE, +, AUU[j,k] * cellvolumes[cell] , dofs_u[k], dofs_u[j])
                    end
                end
                rhs.entries[dofs_u[j]] += bU[j] * cellvolumes[cell]
            end

            fill!(AUU, 0)
            fill!(bU, 0)
        end

        ## apply penalties for boundary condition
        fixed_dofs = boundarydata!(Solution[1], boundary_conditions[1])
        for dof in fixed_dofs
            SE[dof,dof] = fixed_penalty
            rhs.entries[dof] = 0 # fixed_penalty * Solution.entries[dof]
            Solution.entries[dof] = 0
        end

        ## apply global constraints
        for constraint in global_constraints
            apply_constraint!(S, rhs, constraint, Solution)
        end

        flush!(SE)
    end

    if solve_polarisation
        jac_polarisation_view = view(jac_polarisation,indices_∇P,dim^2+1:dim^2+dim)
        APP::Array{Float64,2} = zeros(Float64, ndofs_P, ndofs_P)
        bP::Array{Float64,1} = zeros(Float64, ndofs_P)
    end

    function update_polarisation(; coupled = false)

        if coupled
            offset = FES[1].ndofs
            fill!(S[2,1],0)
            fill!(S[2,2],0)
            fill!(rhs[2],0)
        else
            offset = 0
            fill!(S[1,1],0)
            fill!(rhs[1],0)
        end


        for cell = 1 : ncells

            ## update FE basis evaluators
            item_information[1] = cell
            item_information[2] = cell
            item_information[3] = cellregions[cell]
            pola_op = polarisation_operator4region(item_information)
            ∇u.citem[] = cell
            update_basis!(∇u) 
            ∇P.citem[] = cell
            update_basis!(∇P) 
            dofs_u = view(celldofs_u, : , cell)
            dofs_P = view(celldofs_P, : , cell)

            for qp = 1 : nweights
                ## evaluate ∇u and ∇P in current solution at current quadrature point
                fill!(eval_∇P, 0)
                for j = 1 : ndofs_u, k in indices_∇u
                    eval_∇P[k] += SolutionE[dofs_u[j]] * vals_∇u[k,j,qp] 
                end
                for j = 1 : ndofs_P, k in indices_∇P
                    eval_∇P[9+k] += SolutionE[dofs_P[j]+offsets[2]] * vals_∇P[k,j,qp] 
                end

                ## evaluate jacobian of displacement tensor
                if use_sparsity_detection
                    forwarddiff_color_jacobian!(jac_polarisation, pola_op, eval_∇P, cfg_polarisation)
                else
                    ForwardDiff.vector_mode_jacobian!(Dresult_polarisation,pola_op,result_polarisation,eval_∇P,cfg_polarisation)
                end

                # update matrix
                for j = 1 : ndofs_P
                    # multiply with jacobian
                    mul!(tempP, jac_polarisation_view, view(vals_∇P,:,j,qp))

                    # multiply test function operator evaluation
                    for k = 1 : ndofs_P
                        APP[j,k] += dot(tempP, view(vals_∇P,:,k,qp)) * weights[qp]
                    end
                end 

                # update rhs
                mul!(tempP, jac_polarisation_view, view(eval_∇P, indices_∇P))
                pola_op(value_polarisation, eval_∇P)
                tempP .-= value_polarisation
                for j = 1 : ndofs_P
                    bP[j] += dot(tempP, view(vals_∇P,:,j,qp)) * weights[qp]
                end
            end

            # write into full matrix
            for j = 1 : ndofs_P
                for k = 1 : ndofs_P
                    if abs(APP[j,k]) > 0
                        rawupdateindex!(SE, +, APP[j,k] * cellvolumes[cell] , dofs_P[k], dofs_P[j])
                    end
                end
                rhs.entries[dofs_P[j]] += bP[j] * cellvolumes[cell]
            end

            fill!(APP, 0)
            fill!(bP, 0)
        end

        ## apply penalties for boundary condition
        fixed_dofsP = [offset + 1]
        for dof in fixed_dofsP
            SE[dof,dof] = fixed_penalty
            rhs.entries[dof] = 0 # fixed_penalty * Solution.entries[dof]
            Solution.entries[dof] = 0
        end

        flush!(SE)
    end

    ## prepare other stuff for loop
    emb_params_target = deepcopy(emb_params)
    factorization = nothing
    residual = zeros(Float64, size(SE,1))
    nlres::Float64 = 0
    linres::Float64 = 0
    change::Float64 = 0
    time_assembly::Float64 = 0
    time_solver::Float64 = 0
    iteration::Int = 0

    for step = 1 : nsteps[1]
        emb_params .= nsteps[1] == 1 ? emb_params_target : (step-1)/(nsteps[1]-1) .* emb_params_target

        @info "Entering embedding step = $step/$(nsteps[1]), emb_params = $emb_params"
            
        iteration = 0
        time_solver = 0
        while iteration <= maxiterations[1]
            # update system
            time_assembly = @elapsed update_displacement(; coupled = coupled)
            if iteration == 0 && step == 1
                factorization=linsolver(SE)
                factorize!(factorization, SE)
            end
        
            # compute nonlinear residual
            fill!(residual, 0)
            if coupled
                mul!(residual, SE, Solution.entries)
            else
                mul!(residual, SE, view(Solution[1]))
            end
            residual .-= rhs.entries
            view(residual, fixed_dofs) .= 0

            nlres = norm(residual)
            @info "         ---> nonlinear iteration $iteration, linres = $linres, nlres = $nlres, timeASSEMBLY + timeSOLVE = $time_assembly + $time_solver"
            if nlres < target_residual[1]
                @info "target reached, stopping..."
                break
            end

            # update factorization and solve
            if damping > 0
                SolutionOld.entries .= Solution.entries
            end

            time_solver = @elapsed begin
                update!(factorization)
                ldiv!(view(Solution[1]), factorization, rhs.entries)
            end
            iteration += 1

            fill!(residual, 0)
            if coupled
                mul!(residual, SE, Solution.entries)
            else
                mul!(residual, SE, view(Solution[1]))
            end
            residual .-= rhs.entries
            view(residual, fixed_dofs) .= 0
            linres = norm(residual)

            if damping > 0
                view(Solution[1]) .= (1 - damping) * Solution.entries + damping * SolutionOld.entries
                change = norm(SolutionOld.entries - Solution.entries)
            end

        end
    end

    if !coupled && solve_polarisation
        @info "Entering decoupled polarisation solve..."
        S = FEMatrix(FES[2])       # system matrix
        rhs = FEVector(FES[2])     # system rhs
        SE = S.entries
        residual = zeros(Float64, size(SE,1))

        time_assembly = @elapsed update_polarisation(; coupled = coupled)
        factorization=linsolver(SE)
        factorize!(factorization, SE)

        time_solver = @elapsed begin
            update!(factorization)
            ldiv!(view(Solution[2]), factorization, rhs.entries)
        end
        iteration += 1

        fill!(residual, 0)
        mul!(residual, SE, view(Solution[2]))
        residual .-= rhs.entries
        view(residual, fixed_dofsP) .= 0
        linres = norm(residual)

        @info "         ---> linres = $linres, timeASSEMBLY + timeSOLVE = $time_assembly + $time_solver"
        
    end

    return Solution, nlres

end