### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 57020ab2-b13b-4eb3-aa7a-13176c44c610
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using PlutoUI
    using SimplexGridFactory
	using Triangulate
	using GridVisualize
	using ExtendableGrids
	using PlutoVista
	default_plotter!(PlutoVista)
	TableOfContents()
end

# ╔═╡ 2e3e2ef1-e63a-4ab3-9613-e99f88022d0e
function nanowire_tensorgrid(; scale = [1,1,1,1], nrefs = 1, cut_levels = scale[4]/2, α = nothing, Plotter = nothing, z_levels_dist = 100, version = 1)

    @info "Generating nanowire grid for scale = $scale"

    builder=SimplexGridBuilder(Generator=Triangulate)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ = scale[3]
	if α !== nothing
	    vol_factor_core = 4.0^-1
    	vol_factor_shell = 4.0^-1
	else
		vol_factor_core = 4.0^-(nrefs-1)
    	vol_factor_shell = 4.0^-nrefs
	end
    vol_factor_stressor = 4.0^-nrefs
    hz_factor = 2.0^-nrefs

    A_core = 3*sqrt(3)/2 * scale[1]^2
    A_shell = 3*sqrt(3)/2 * (scale[2]^2 + 2*scale[1]*scale[2])
    A_stressor = sqrt(3)/2 * scale[3] * (7*(scale[1]+scale[2]) + 3*scale[3])
    if α !== nothing
        A_interface = 3*(d2 * sqrt(3)/2*α)
        A_shell = A_shell - A_interface
        A_stressor = A_stressor - A_interface
        vol_factor_interface = 4.0^-nrefs
    end

    # bottom side at Z = 0
	#p99= point!(builder,0,-d2)
    p0 = point!(builder,0,0)
    p1 = point!(builder,d1,0)
    p2 = point!(builder,d1/2,sqrt(3)/2*d1)
    p3 = point!(builder,-d1/2,sqrt(3)/2*d1)
    p4 = point!(builder,-d1,0)
    p5 = point!(builder,-d1/2,-sqrt(3)/2*d1)
    p6 = point!(builder,d1/2,-sqrt(3)/2*d1)
    p7 = point!(builder,d2,0)
    p8 = point!(builder,d2/2,sqrt(3)/2*d2)
    p9 = point!(builder,-d2/2,sqrt(3)/2*d2)
    p10 = point!(builder,-d2,0)
    p11 = point!(builder,-d2/2,-sqrt(3)/2*d2)
    p12 = point!(builder,d2/2,-sqrt(3)/2*d2)
    if version == 1
        p13 = point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ)
        p14 = point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ)
    elseif version == 2
        p13 = point!(builder,d2,-δ)
        p14 = point!(builder,d2/2,-sqrt(3)/2*d2-δ)
        p15 = point!(builder,-d2/2,-sqrt(3)/2*d2-δ)
        p16 = point!(builder,-d2,-δ)
    end

    if α !== nothing
        xstar = α/2 + d2*(1 - sqrt(3)*α/(4*δ))
        ystar = - α*(2*sqrt(3)*δ + 3*d2)/(4*δ)
        p15 = point!(builder,-d2+α/2,sqrt(3)/2*α)
        p16 = point!(builder,(α-d2)/2,-sqrt(3)/2*(d2-α))
        p17 = point!(builder,(d2-α)/2,-sqrt(3)/2*(d2-α))
        p18 = point!(builder,d2-α/2,sqrt(3)/2*α)
        p19 = point!(builder,xstar,ystar)
        p20 = point!(builder,(d2+α)/2,-sqrt(3)/2*(d2+α))
        p21 = point!(builder,-(d2+α)/2,-sqrt(3)/2*(d2+α))
        p22 = point!(builder,-xstar,ystar)
    end

    if α !== nothing
        facetregion!(builder,1) # core region
        facet!(builder,p1,p2)
        facet!(builder,p0,p1)
        facet!(builder,p0,p2)

        facet!(builder,p2,p3)
        facet!(builder,p0,p3)

        facet!(builder,p3,p4)
        facet!(builder,p0,p4)

        facet!(builder,p4,p5)
        facet!(builder,p0,p5)

        facet!(builder,p5,p6)
        facet!(builder,p0,p6)

        facet!(builder,p6,p1)

        facetregion!(builder,2) # shell region
        facet!(builder,p15,p16)
        facet!(builder,p16,p17)
        facet!(builder,p17,p18)
        facet!(builder,p18,p8)
        facet!(builder,p8,p9)
        facet!(builder,p9,p15)

        facetregion!(builder,3) # interface region (inside)
        facet!(builder,p15,p10)
        facet!(builder,p10,p11)
        facet!(builder,p11,p12)
        facet!(builder,p12,p7)
        facet!(builder,p7,p18)

        facetregion!(builder,4) # interface region (outside)
        facet!(builder,p7,p19)
        facet!(builder,p19,p20)
        facet!(builder,p20,p21)
        facet!(builder,p21,p22)
        facet!(builder,p22,p10)

        facetregion!(builder,5) # stressor region
        facet!(builder,p19,p13)
        facet!(builder,p13,p14)
        facet!(builder,p14,p22)

        cellregion!(builder,1) # material 1 (core)
        maxvolume!(builder,A_core/12*vol_factor_core)
        regionpoint!(builder,(d1/2,sqrt(3)/4*d1))
        regionpoint!(builder,(0,sqrt(3)/4*d1))
        regionpoint!(builder,(-d1/2,sqrt(3)/4*d1))
        regionpoint!(builder,(d1/2,-sqrt(3)/4*d1))
        regionpoint!(builder,(0,-sqrt(3)/4*d1))
        regionpoint!(builder,(-d1/2,-sqrt(3)/4*d1))

        cellregion!(builder,2) # material 2 (shell)
        maxvolume!(builder,A_shell*vol_factor_shell)
        regionpoint!(builder,(scale[1]+scale[2]/2,0))

        cellregion!(builder,3) # material 3 (inside interface)
        maxvolume!(builder,A_interface*vol_factor_interface/2)
        regionpoint!(builder,(d2-α/2,0))

        cellregion!(builder,4) # material 4 (outside interface)
        maxvolume!(builder,A_interface*vol_factor_interface/2)
        regionpoint!(builder,(0,-sqrt(3)/2*(d2+α/2)))

        cellregion!(builder,5) # material 5 (stressor)
        maxvolume!(builder,A_stressor*vol_factor_stressor)
        regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2))
    else
        facetregion!(builder,1) # core region
        facet!(builder,p1,p2)
        facet!(builder,p2,p3)
        facet!(builder,p3,p4)
        facet!(builder,p4,p5)
        facet!(builder,p5,p6)
        facet!(builder,p6,p1)

        facetregion!(builder,2) # shell region
        facet!(builder,p7,p8)
        facet!(builder,p8,p9)
        facet!(builder,p9,p10)
        facet!(builder,p10,p11)
        facet!(builder,p11,p12)
        facet!(builder,p12,p7)

        facetregion!(builder,3) # stressor region
        if version == 1
            facet!(builder,p7,p13)
            facet!(builder,p13,p14)
            facet!(builder,p14,p10)
        elseif version == 2
            facet!(builder,p7,p13)
            facet!(builder,p13,p14)
            facet!(builder,p14,p15)
            facet!(builder,p15,p16)
            facet!(builder,p16,p10)
        end

        cellregion!(builder,1) # material 1
        maxvolume!(builder,A_core/6*vol_factor_core)
        regionpoint!(builder,(0,0))

        cellregion!(builder,2) # material 2
        maxvolume!(builder,A_shell/6*vol_factor_shell)
        regionpoint!(builder,(scale[1]+scale[2]/2,0))

        cellregion!(builder,3) # material 3
        maxvolume!(builder,A_stressor/6*vol_factor_stressor)
        regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2))
    end
    xgrid = simplexgrid(builder)

    hz = z_levels_dist * hz_factor
    if cut_levels !== nothing

        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
        for i = 1 : length(cut_levels)
            index = findfirst(item -> item >= cut_levels[i], z_levels)
            if z_levels_nonuniform[index] == cut_levels[i]
                z_levels_nonuniform[index] = cut_levels[i]-hz/2:hz/4:cut_levels[i]+hz/2
            else
                hz1 = (cut_levels[i]-z_levels_nonuniform[index-1])/2
                hz2 = (z_levels_nonuniform[index]-cut_levels[i])/2
                z_levels_nonuniform[index-1] = z_levels_nonuniform[index-1]:hz1:cut_levels[i]
                z_levels_nonuniform[index]   = cut_levels[i]+hz2:hz2:z_levels_nonuniform[index]
            end
        end
    else
        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
    end
    z_levels_nonuniform = vcat(z_levels_nonuniform...)

    if α !== nothing
        cellregions = xgrid[CellRegions]
        for i = 1 : num_cells(xgrid)
            if cellregions[i] == 3
                cellregions[i] = 2
            end
            if cellregions[i] == 4 || cellregions[i] == 5
                cellregions[i] = 3
            end
        end
        xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset = 5, top_offset = 8)
        # the offsets lead to the following boundary regions:
        # 1 = side core (not seen from outside)
        # 2 = side shell
        # 3 = side interface inside
        # 4 = side interface outside
        # 5 = side stressor
        # 6 = bottom core
        # 7 = bottom shell
        # 8 = bottom stressor
        # 9 = top core
        # 10 = top shell
        # 11 = top stressor
    else
        xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset = 3, top_offset = 6)
        # the offsets lead to the following boundary regions:
        # 1 = side core (not seen from outside)
        # 2 = side shell
        # 3 = side stressor
        # 4 = bottom core
        # 5 = bottom shell
        # 6 = bottom stressor
        # 7 = top core
        # 8 = top shell
        # 9 = top stressor
    end

    return xgrid
end

# ╔═╡ 6da44973-6c89-4cbc-b15e-33ea72037613
geometry = [30, 20, 10, 2000]

# ╔═╡ 0ac0298e-d786-4623-a89b-3b466f3786d8
xgrid = nanowire_tensorgrid(; scale = geometry, cut_levels = geometry[4]/2, nrefs = 2, α = geometry[3]/4, version = 1)

# ╔═╡ 38b899e7-ca56-41ec-b195-ceaad78072df
	gridplot(xgrid; Plotter=PlutoVista)

# ╔═╡ Cell order:
# ╠═57020ab2-b13b-4eb3-aa7a-13176c44c610
# ╠═2e3e2ef1-e63a-4ab3-9613-e99f88022d0e
# ╠═6da44973-6c89-4cbc-b15e-33ea72037613
# ╠═0ac0298e-d786-4623-a89b-3b466f3786d8
# ╠═38b899e7-ca56-41ec-b195-ceaad78072df
