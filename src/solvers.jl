
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
    Solution = FEVector{Float64}(Problem.unknown_names, FES)

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
            residual = solve!(Solution, Problem; subiterations = subiterations[s], show_statistics = true, linsolver = linsolver, maxiterations = maxiterations[s], target_residual = target_residual[s])
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
            damping_step = 0.05,                                 # checks damping values within the range 0:damping_step:1
            maxiterations = 100) where {Tv,Ti}                  # maximal iterations

    ## generate solution vector
    FES = Array{FESpace{Tv,Ti},1}(undef, length(FETypes))
    for j = 1 : length(FETypes)
        FES[j] = FESpace{FETypes[j]}(xgrid)
    end
    Solution = FEVector{Float64}(Problem.unknown_names, FES)
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