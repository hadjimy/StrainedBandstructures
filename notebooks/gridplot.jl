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
	default_plotter!(PlutoVista)
	TableOfContents()
end

# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole NanoWiresJulia project."

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

# ╔═╡ baf477b5-a101-4a03-83a9-860540398a6e
md"""
## bimetal\_tensorgrid\_uniform
"""

# ╔═╡ f6c827e6-b004-4c9b-b560-ddc609631b64
md"""
# Nanowire 3D
"""

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

# ╔═╡ 4cdc651c-54ea-4645-b8e1-0202c47e7f2c
grid_bimetal_tensorgrid=bimetal_tensorgrid()

# ╔═╡ 0af4926e-eb19-4dca-9773-ecbc1e109c6a
md"""
## nanowire_tensorgrid
Three-region nanowire generated with 2D triangulation and tensor expansion in the z-direction.
"""

# ╔═╡ 9d967453-a5ff-4fb9-8346-5a1e734e59b2
begin
	geometry = [30, 20, 10, 500]
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
end

# ╔═╡ ee5a4369-8ed2-44fd-b4f1-2a19f485abfb
begin
	grid_nanowire_tensorgrid1 = nanowire_tensorgrid(; scale = scale, cut_levels = geometry[4]/2, nrefs = 1, α = geometry[3]/4, version = 1)
	gridplot(grid_nanowire_tensorgrid1)
end

# ╔═╡ b93e852a-ed8e-4250-a240-aca1da136184
begin
	grid_nanowire_tensorgrid2 = nanowire_tensorgrid(; scale = scale, cut_levels = geometry[4]/2, nrefs = 2, α = geometry[3]/2, version = 1)
	gridplot(grid_nanowire_tensorgrid2)
end

# ╔═╡ cf0ca13b-b57b-4263-adf9-6d168ae75f0d
gridplot(grid_bimetal_tensorgrid)

# ╔═╡ 7dcb4899-5199-4201-82ee-64b046ea68ce
grid_bimetal_tensorgrid_uniform=bimetal_tensorgrid_uniform(nrefs=2)

# ╔═╡ Cell order:
# ╟─882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
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
# ╟─baf477b5-a101-4a03-83a9-860540398a6e
# ╟─f6c827e6-b004-4c9b-b560-ddc609631b64
# ╟─1c4624f1-d56c-47c5-869b-c29a4aab3f06
# ╠═ece590ef-61ad-47b6-9888-4e4c0d8163e7
# ╠═4cdc651c-54ea-4645-b8e1-0202c47e7f2c
# ╟─0af4926e-eb19-4dca-9773-ecbc1e109c6a
# ╠═9d967453-a5ff-4fb9-8346-5a1e734e59b2
# ╠═ee5a4369-8ed2-44fd-b4f1-2a19f485abfb
# ╠═b93e852a-ed8e-4250-a240-aca1da136184
# ╠═cf0ca13b-b57b-4263-adf9-6d168ae75f0d
# ╠═7dcb4899-5199-4201-82ee-64b046ea68ce
