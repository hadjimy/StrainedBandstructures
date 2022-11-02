### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using PlutoUI
	using NanoWiresJulia
	using ExtendableGrids
	using GridVisualize
	using PlutoVista
	using Triangulate
	using SimplexGridFactory
	using LinearAlgebra
	default_plotter!(PlutoVista)
	TableOfContents()
end

# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole NanoWiresJulia project."

# ╔═╡ f6c827e6-b004-4c9b-b560-ddc609631b64
md"""
# Nanowire 3D
"""

# ╔═╡ 0af4926e-eb19-4dca-9773-ecbc1e109c6a
md"""
## nanowire_tensorgrid
Three-region nanowire generated with 2D triangulation and tensor expansion in the z-direction.
"""

# ╔═╡ ba81b798-9a56-450e-9f64-0f93f5c01913
function nanowire_tensorgrid!(; scale = [1,1,1,1], nrefs = 1,
	cut_levels = scale[4]/2, α = nothing, Plotter = nothing, z_levels_dist = 100, version = 1, corner_refinement = false, manual_refinement = false)
	
    @info "Generating nanowire grid for scale = $scale"

    builder = SimplexGridBuilder(Generator=Triangulate)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ  = scale[3]

    A_core = 3*sqrt(3)/2 * scale[1]^2
    A_shell = 3*sqrt(3)/2 * (scale[2]^2 + 2*scale[1]*scale[2])
    A_stressor = sqrt(3)/2 * scale[3] * (7*(scale[1]+scale[2]) + 3*scale[3])
    if α !== nothing
        A_interface = 3*(d2 * sqrt(3)/2*α)
        A_shell = A_shell - A_interface
        A_stressor = A_stressor - A_interface

        vol_factor_core = 4.0^0
        vol_factor_shell = 4.0^0
        vol_factor_interface = 4.0^-(nrefs+1)
        vol_factor_stressor = 4.0^0
    else
        vol_factor_core = 4.0^-nrefs
        vol_factor_shell = 4.0^-nrefs
        vol_factor_stressor = 4.0^-nrefs
    end
    hz_factor = 2.0^-nrefs

    # bottom side at Z = 0
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
		if manual_refinement == true
			num_pts = 8
			xstar = α/4 + d2*(1 - sqrt(3)*α/2/(4*δ))
        	ystar = - α/2*(2*sqrt(3)*δ + 3*d2)/(4*δ)
			for n = 0 : num_pts
				g = n/num_pts
				## adding points along the horizontal interface
				# convex combination between p11 & p16 midpoint and p12 & p17 midpoint
				px = g*(α/2-d2)/2 + (1-g)*(d2-α/2)/2
				py = -sqrt(3)/2*(d2-α/2)
				point!(builder,px,py)
				# convex combination between p11 & p21 midpoint and p12 & p20 midpoint
				px = g*(-(d2+α/2)/2) + (1-g)*(d2+α/2)/2
				py = -sqrt(3)/2*(d2+α/2)
				point!(builder,px,py)

				## adding points along the left interface
				# convex combination between p11 & p16 midpoint and p10 & p15 midpoint
				px = g*(α/2-d2)/2 + (1-g)*(-d2+α/4)
				py = g*(-sqrt(3)/2*(d2-α/2)) + (1-g)*sqrt(3)/2*α/2
				point!(builder,px,py)
				# convex combination between p11 & p21 midpoint and p10 & p22 midpoint
				px = g*(-(d2+α/2)/2) + (1-g)*(-xstar)
				py = g*(-sqrt(3)/2*(d2+α/2)) + (1-g)*ystar
				point!(builder,px,py)

				## adding points along the right interface
				# convex combination between p12 & p17 midpoint and p7 & p18 midpoint
				px = g*(d2-α/2)/2 + (1-g)*(d2-α/4)
				py = g*(-sqrt(3)/2*(d2-α/2)) + (1-g)*(sqrt(3)/2*α/2)
				point!(builder,px,py)
				# convex combination between p12 & p20 midpoint and p7 & p19 midpoint
				px = g*(d2+α/2)/2 + (1-g)*(xstar)
				py = g*(-sqrt(3)/2*(d2+α/2))+(1-g)*ystar
				point!(builder,px,py)
			end
		end
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
        maxvolume!(builder,A_interface*vol_factor_interface)
        regionpoint!(builder,(d2-α/2,0))

        cellregion!(builder,4) # material 4 (outside interface)
        maxvolume!(builder,A_interface*vol_factor_interface)
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

	if corner_refinement == true
	    function unsuitable(x1,y1,x2,y2,x3,y3, area)
	        bary = [(x1+x2+x3)/3,(y2+y2+y3)/3]
	        dist = norm(bary-refinement_center)
	        if area > 1.0*dist
	            return 1
	        else
	            return 0
	        end
	    end
	
		refinement_center = [-d2/2,-sqrt(3)/2*d2]
	    options!(builder, unsuitable=unsuitable)
		xgrid1 = simplexgrid(builder)
	
		# refinement_center = [d2/2,-sqrt(3)/2*d2]
		# options!(builder, unsuitable=unsuitable)
		# xgrid2 = simplexgrid(builder)
		
		# xgrid = glue(xgrid1,xgrid2)
		xgrid = xgrid1
	else
		xgrid = simplexgrid(builder)
	end
	
    hz = z_levels_dist * hz_factor
    if cut_levels !== nothing

        z_levels = 0:hz:scale[4]
        z_levels_nonuniform = Vector{Any}(z_levels)
        for i = 1 : length(cut_levels)
            index = findfirst(item -> item >= cut_levels[i], z_levels)
            if z_levels_nonuniform[index] == cut_levels[i]
                z_levels_nonuniform[index] =
					cut_levels[i]-hz/2:hz/4:cut_levels[i]+hz/2
            else
                hz1 = (cut_levels[i]-z_levels_nonuniform[index-1])/2
                hz2 = (z_levels_nonuniform[index]-cut_levels[i])/2
                z_levels_nonuniform[index-1] =
					z_levels_nonuniform[index-1]:hz1:cut_levels[i]
                z_levels_nonuniform[index] =
					cut_levels[i]+hz2:hz2:z_levels_nonuniform[index]
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
        xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=5, top_offset=8)
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
        xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=3, top_offset=6)
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

# ╔═╡ 9d967453-a5ff-4fb9-8346-5a1e734e59b2
begin
	geometry = [30, 20, 10, 500]
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
end

# ╔═╡ ee5a4369-8ed2-44fd-b4f1-2a19f485abfb
begin
	grid_nanowire_tensorgrid1 = nanowire_tensorgrid!(; scale = scale,
		cut_levels = geometry[4]/2, nrefs = 2, α = geometry[3]/2, version = 1,
		corner_refinement = false, manual_refinement = true)
	gridplot(grid_nanowire_tensorgrid1)
end

# ╔═╡ 52dd26ec-5f74-4214-b3f9-5b0c84277099
begin
	grid_nanowire_tensorgrid_corners = nanowire_tensorgrid!(; scale = scale, cut_levels = geometry[4]/2, nrefs = 2, α = geometry[3]/2, version = 1, corner_refinement = true, manual_refinement = true)
	gridplot(grid_nanowire_tensorgrid_corners)
end

# ╔═╡ b93e852a-ed8e-4250-a240-aca1da136184
begin
	grid_nanowire_tensorgrid2 = nanowire_tensorgrid(; scale = scale, cut_levels = geometry[4]/2, nrefs = 1, α = geometry[3]/4, version = 1)
	gridplot(grid_nanowire_tensorgrid2)
end

# ╔═╡ 40b2766e-d66e-4e07-9982-803392c35dee
begin
	grid_nanowire_tensorgrid3 = nanowire_tensorgrid(; scale = scale,
		cut_levels = geometry[4]/2, nrefs = 2, α = geometry[3]/2, version = 1)
	gridplot(grid_nanowire_tensorgrid3)
end

# ╔═╡ 1c4624f1-d56c-47c5-869b-c29a4aab3f06
md"""
## nanowire_grid
Three-region nanowire generated with TetGen
"""

# ╔═╡ ece590ef-61ad-47b6-9888-4e4c0d8163e7
begin
	grid_nanowire_grid = nanowire_grid(reflevel=1)
	gridplot(grid_nanowire_grid)
end

# ╔═╡ 258594c5-b5fc-4720-ab59-3e0e71826df0
md"""
# Bimetal 2D
"""

# ╔═╡ e1d32a8e-49be-4c8e-9218-c01843368f73
md"""
## bimetal_strip2D
"""

# ╔═╡ 2608784c-d93e-41f6-a2b4-b2fcd63d8dd8
begin
	grid_bimetal_strip2D=bimetal_strip2D(scale=[50,100])
	gridplot(grid_bimetal_strip2D)
end

# ╔═╡ 44a96db2-0faa-4536-a062-317d7e967f4f
md"""
## condensator2D
"""

# ╔═╡ a786e949-aa71-475d-a82c-9190936b0f51
begin
	grid_condensator2D=condensator2D(A = 50, B = 100, d = 5, reflevel = 1)
	gridplot(grid_condensator2D)
end

# ╔═╡ e1ac16e9-a174-4a2d-8a66-a1e8eb807732
md"""
# Bimetal 3D
"""

# ╔═╡ b2d93b2f-78c9-4a4d-8cac-b04a625a3f75
md"""
## bimetal_strip3D
"""

# ╔═╡ 0b76e405-efd4-4c0b-8d53-486cf164e04c
grid_bimetal_strip3D=bimetal_strip3D(reflevel=2)

# ╔═╡ ca3fb071-ffda-46fc-ac22-0f3ecabca06a
gridplot(grid_bimetal_strip3D,xplanes=0.5,yplanes=0.5,zplanes=0.5)

# ╔═╡ ef096445-8cf2-4d48-8011-e369527ad103
md"""
## bimetal\_strip3D\_middle_layer
"""

# ╔═╡ d1c2cc0a-2d09-46b2-88b1-4c8358f42d2a
grid_bimetal_strip3D_middle_layer=bimetal_strip3D_middle_layer()

# ╔═╡ 5048badf-f655-4e8c-9455-c14a6913c53b
unique(grid_bimetal_strip2D[BFaceRegions])

# ╔═╡ 3f648b46-dd2c-4ffd-9374-2b2493ac3d0c
md"""
## condensator3D_tensorgrid
"""

# ╔═╡ 7c1178fc-b6db-48dd-83c0-cf8a266f2655
grid_condensator3D_tensorgrid=condensator3D_tensorgrid()

# ╔═╡ 7c54868c-756c-4a52-a07a-451d0677f6c1
gridplot(grid_condensator3D_tensorgrid,xplanes=[25.0])

# ╔═╡ f78a3b26-a222-44f0-96b7-529e71dc0e81
md"""
## bimetal_tensorgrid
"""

# ╔═╡ 4cdc651c-54ea-4645-b8e1-0202c47e7f2c
grid_bimetal_tensorgrid=bimetal_tensorgrid()

# ╔═╡ baf477b5-a101-4a03-83a9-860540398a6e
md"""
## bimetal\_tensorgrid\_uniform
"""

# ╔═╡ cf0ca13b-b57b-4263-adf9-6d168ae75f0d
gridplot(grid_bimetal_tensorgrid)

# ╔═╡ 7dcb4899-5199-4201-82ee-64b046ea68ce
grid_bimetal_tensorgrid_uniform=bimetal_tensorgrid_uniform(nrefs=2)

# ╔═╡ 8bef89d0-a429-4571-9bc4-e4aef9a0a536
function square_localref()

    builder=SimplexGridBuilder(Generator=Triangulate)
    cellregion!(builder,1)
    maxvolume!(builder,0.01)
    regionpoint!(builder,0.5,0.5)

    p1=point!(builder,0,0)
    p2=point!(builder,1,0)
    p3=point!(builder,1,1)
    p4=point!(builder,0,1)

    facetregion!(builder,1)
    facet!(builder,p1,p2)
    facetregion!(builder,2)
    facet!(builder,p2,p3)
    facetregion!(builder,3)
    facet!(builder,p3,p4)
    facetregion!(builder,4)
    facet!(builder,p4,p1)

    refinement_center=[0.5,0.5]
    function unsuitable(x1,y1,x2,y2,x3,y3, area)
        bary=[(x1+x2+x3)/3,(y2+y2+y3)/3]
        dist=norm(bary-refinement_center)
        if area > 0.1*dist
            return 1
        else
            return 0
        end
    end
    options!(builder, unsuitable=unsuitable)

	xgrid = simplexgrid(builder)
	
	return xgrid
end

# ╔═╡ b889dec0-b75b-4188-8b31-c5043228da4a
begin
	grid_square_localref=square_localref()
	gridplot(grid_square_localref)
end

# ╔═╡ Cell order:
# ╟─882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─f6c827e6-b004-4c9b-b560-ddc609631b64
# ╟─0af4926e-eb19-4dca-9773-ecbc1e109c6a
# ╟─ba81b798-9a56-450e-9f64-0f93f5c01913
# ╠═9d967453-a5ff-4fb9-8346-5a1e734e59b2
# ╠═ee5a4369-8ed2-44fd-b4f1-2a19f485abfb
# ╠═52dd26ec-5f74-4214-b3f9-5b0c84277099
# ╠═b93e852a-ed8e-4250-a240-aca1da136184
# ╠═40b2766e-d66e-4e07-9982-803392c35dee
# ╟─1c4624f1-d56c-47c5-869b-c29a4aab3f06
# ╠═ece590ef-61ad-47b6-9888-4e4c0d8163e7
# ╟─258594c5-b5fc-4720-ab59-3e0e71826df0
# ╟─e1d32a8e-49be-4c8e-9218-c01843368f73
# ╠═2608784c-d93e-41f6-a2b4-b2fcd63d8dd8
# ╟─44a96db2-0faa-4536-a062-317d7e967f4f
# ╟─a786e949-aa71-475d-a82c-9190936b0f51
# ╟─e1ac16e9-a174-4a2d-8a66-a1e8eb807732
# ╟─b2d93b2f-78c9-4a4d-8cac-b04a625a3f75
# ╠═0b76e405-efd4-4c0b-8d53-486cf164e04c
# ╠═ca3fb071-ffda-46fc-ac22-0f3ecabca06a
# ╟─ef096445-8cf2-4d48-8011-e369527ad103
# ╠═d1c2cc0a-2d09-46b2-88b1-4c8358f42d2a
# ╠═5048badf-f655-4e8c-9455-c14a6913c53b
# ╟─3f648b46-dd2c-4ffd-9374-2b2493ac3d0c
# ╠═7c1178fc-b6db-48dd-83c0-cf8a266f2655
# ╠═7c54868c-756c-4a52-a07a-451d0677f6c1
# ╟─f78a3b26-a222-44f0-96b7-529e71dc0e81
# ╠═4cdc651c-54ea-4645-b8e1-0202c47e7f2c
# ╟─baf477b5-a101-4a03-83a9-860540398a6e
# ╠═cf0ca13b-b57b-4263-adf9-6d168ae75f0d
# ╠═7dcb4899-5199-4201-82ee-64b046ea68ce
# ╠═8bef89d0-a429-4571-9bc4-e4aef9a0a536
# ╠═b889dec0-b75b-4188-8b31-c5043228da4a
