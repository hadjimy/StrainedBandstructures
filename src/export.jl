
## writes (upscaled) Displacement and (optionally) Polarisation into vtk file
## on unbended and bended grid
function exportVTK(filename, Displacement::FEVectorBlock{T,Tv,Ti,FEType,APT}, Polarisation = nothing;
    eps_gfind = 1e-12,
    upscaling = 0,
    P0strain::Bool = true,
    strain_model = NonlinearStrain2D) where {T,Tv,Ti,FEType,APT}

    ## get original grid
    xgrid_plot = Displacement.FES.xgrid
    ncomponents = get_ncomponents(FEType)

    ## prepare strain interpolation
    size_ϵu = Int((ncomponents+1)*ncomponents/2)
    function interpolate_postprocess(result, input)
        # input = [∇u]
        eval_strain!(result, input, strain_model)
        return nothing
    end
    interpolate_action = Action(interpolate_postprocess, [size_ϵu, ncomponents^2])

    ## upscaling (= interpolate to P1 on finer grid)
    Solution_plot = nothing
    if upscaling > 0
        xgrid_plot = uniform_refine(xgrid_plot, upscaling; store_parents = true)
        FES = FESpace{H1P1{ncomponents}}(xgrid_plot)
        FES_strain = FESpace{P0strain ? L2P0{size_ϵu} : H1P1{size_ϵu}}(xgrid_plot)
        if Polarisation !== nothing
            FES_polarisation = FESpace{H1P1{1}}(xgrid_plot)
            Solution_plot = FEVector{Float64}(["u_h (upscale)", "ϵu", "V_P"],[FES, FES_strain, FES_polarisation])
            interpolate!(Solution_plot[3], Polarisation; use_cellparents = true, eps = eps_gfind)
        else
            Solution_plot = FEVector{Float64}(["u_h (upscale)", "ϵu"],[FES, FES_strain])
        end
        interpolate!(Solution_plot[1], Displacement; use_cellparents = true, eps = eps_gfind)
        interpolate!(Solution_plot[2], Displacement; operator = Gradient, postprocess = interpolate_action, use_cellparents = true, eps = eps_gfind)
    else
        Strain = FEVector("ϵu",FESpace{P0strain ? L2P0{size_ϵu} : H1P1{size_ϵu}}(xgrid_plot))
        interpolate!(Strain[1], Displacement; operator = Gradient, postprocess = interpolate_action, eps = eps_gfind)
        if Polarisation !== nothing
            Solution_plot = [Displacement, Strain[1], Polarisation]
        else
            Solution_plot = [Displacement, Strain[1]]
        end
    end

    ## write solution to vtk (unbended and bended)
    kwargs = Dict()
    kwargs[:cellregions] = xgrid_plot[CellRegions]
    #views_displacements = nodevalues_view(Solution_plot[1])
    #kwargs[:displacement_X] = views_displacements[1]
    #kwargs[:displacement_Y] = views_displacements[2]
    #kwargs[:displacement_Z] = views_displacements[3]
    kwargs[:displacement] = view(nodevalues(Solution_plot[1], Identity),:,:)
    kwargs[:grad_displacement] = view(nodevalues(Solution_plot[1], Gradient),:,:)
    kwargs[:strain] = view(nodevalues(Solution_plot[2], Identity),:,:)
    if Polarisation !== nothing
        kwargs[:polarisation] = nodevalues_view(Solution_plot[3])[1]
    end
    ExtendableGrids.writeVTK(filename * "_unbend.vtu", xgrid_plot; kwargs...)
    xgrid_displaced = displace_mesh(xgrid_plot, Solution_plot[1])
    ExtendableGrids.writeVTK(filename * "_bend.vtu", xgrid_displaced; kwargs...)
end