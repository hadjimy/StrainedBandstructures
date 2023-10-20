
## writes (upscaled) Displacement and (optionally) Polarisation into vtk file
## on unbended and bended grid
function exportVTK(filename, Displacement::FEVectorBlock{T,Tv,Ti,FEType,APT}, Polarisation = nothing;
    eps_gfind = 1e-12,
    eps0 = nothing,         # misfit strain (to be subtracted from strain to get elastic strain)
    upscaling = 0,
    P0strain::Bool = true,
    strain_model = NonlinearStrain2D) where {T,Tv,Ti,FEType,APT}

    ## get original grid
    xgrid_plot = Displacement.FES.xgrid
    ncomponents = get_ncomponents(FEType)

    ## prepare strain interpolation
    size_ϵu = Int((ncomponents+1)*ncomponents/2)
    misfit = [[0,0,0],[0,0,0],[0,0,0]]
    function interpolate_postprocess(result, input, qpinfo)
        region = qpinfo.region
        # input = [∇u]
        eval_strain!(result, input, strain_model)
        if typeof(misfit[region]) <: Real
            result[1] -= misfit[region]
            result[2] -= misfit[region]
            result[3] -= misfit[region]
        else
            result[1] -= misfit[region][1]
            result[2] -= misfit[region][2]
            result[3] -= misfit[region][3]
        end
        # result cointains strain in Voigt notation; divide by 2 to calculate the off-diagonal strain components
        result[4] /= 2
        result[5] /= 2
        result[6] /= 2
        return nothing
    end

    ## upscaling (= interpolate to P1 on finer grid)
    Solution_plot = nothing
    if upscaling > 0
        xgrid_plot = uniform_refine(xgrid_plot, upscaling; store_parents = true)
        FES = FESpace{H1P1{ncomponents}}(xgrid_plot)
        FES_strain = FESpace{P0strain ? L2P0{size_ϵu} : H1P1{size_ϵu}}(xgrid_plot)
        if Polarisation !== nothing
            FES_polarisation = FESpace{H1P1{1}}(xgrid_plot)
            Solution_plot = FEVector([FES, FES_strain, FES_strain, FES_polarisation])
            interpolate!(Solution_plot[4], Polarisation; use_cellparents = true, eps = eps_gfind)
        else
            Solution_plot = FEVector([FES, FES_strain, FES_strain])
        end
        lazy_interpolate!(Solution_plot[1], [Displacement], [id(1)]; use_cellparents = true, eps = eps_gfind)
        lazy_interpolate!(Solution_plot[2], [Displacement], [grad(1)]; postprocess = interpolate_postprocess, use_cellparents = true, eps = eps_gfind)
        if eps0 !== nothing
            misfit = eps0
            lazy_interpolate!(Solution_plot[3], [Displacement], [grad(1)]; postprocess = interpolate_postprocess, use_cellparents = true, eps = eps_gfind)
        else
            lazy_interpolate!(Solution_plot[3], [Displacement], [grad(1)]; postprocess = interpolate_postprocess, use_cellparents = true, eps = eps_gfind)
        end
    else
        Strain = FEVector(FESpace{P0strain ? L2P0{size_ϵu} : H1P1{size_ϵu}}(xgrid_plot))
        lazy_interpolate!(Strain[1], [Displacement], [grad(1)]; postprocess = interpolate_postprocess, eps = eps_gfind)
        if eps0 !== nothing
            Strain2 = FEVector(FESpace{P0strain ? L2P0{size_ϵu} : H1P1{size_ϵu}}(xgrid_plot))
            misfit = eps0
            lazy_interpolate!(Strain2[1], [Displacement], [grad(1)]; postprocess = interpolate_postprocess, eps = eps_gfind)
        else
            Strain2 = Strain
        end
        if Polarisation !== nothing
            Solution_plot = [Displacement, Strain[1], Strain2[1], Polarisation]
        else
            Solution_plot = [Displacement, Strain[1], Strain2[1]]
        end
    end

    ## write solution to vtk (unbended and bended)
    kwargs = Dict()
    kwargs[:cellregions] = xgrid_plot[CellRegions]
    kwargs[:displacement] = view(nodevalues(Solution_plot[1], Identity),:,:)
    kwargs[:grad_displacement] = view(nodevalues(Solution_plot[1], Gradient),:,:)
    kwargs[:strain] = view(nodevalues(Solution_plot[2], Identity),:,:)
    if eps0 !== nothing
        kwargs[:elastic_strain] = view(nodevalues(Solution_plot[3], Identity),:,:)
    end
    if Polarisation !== nothing
        kwargs[:polarisation] = nodevalues_view(Solution_plot[4])[1]
    end
    ExtendableGrids.writeVTK(filename * "_unbend.vtu", xgrid_plot; kwargs...)
    xgrid_displaced = displace_mesh(xgrid_plot, Solution_plot[1])
    ExtendableGrids.writeVTK(filename * "_bend.vtu", xgrid_displaced; kwargs...)
end