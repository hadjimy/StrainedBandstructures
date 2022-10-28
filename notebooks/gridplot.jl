### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
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

# ╔═╡ b8ec1d59-8a87-4cf0-bde8-790b420f71bc


# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole NanoWiresJulia project."

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

# ╔═╡ e1d32a8e-49be-4c8e-9218-c01843368f73
md"""
## bimetal_strip2D
"""

# ╔═╡ 2608784c-d93e-41f6-a2b4-b2fcd63d8dd8
grid_bimetal_strip2D=bimetal_strip2D()

# ╔═╡ 2ede7b7b-1b03-49df-99fb-9137a877ab22
gridplot(grid_bimetal_strip2D)

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

# ╔═╡ 44a96db2-0faa-4536-a062-317d7e967f4f
md"""
## condensator2D
"""

# ╔═╡ a786e949-aa71-475d-a82c-9190936b0f51
grid_condensator2D=condensator2D()

# ╔═╡ 719476ce-60c6-42ec-bf02-4679a0e67edc
gridplot(grid_condensator2D)

# ╔═╡ 1c4624f1-d56c-47c5-869b-c29a4aab3f06
md"""
## nanowire_grid
"""

# ╔═╡ ece590ef-61ad-47b6-9888-4e4c0d8163e7
grid_nanowire_grid=nanowire_grid()

# ╔═╡ f85c2e83-d151-4647-a88e-15225c62631d
gridplot(grid_nanowire_grid)

# ╔═╡ 0af4926e-eb19-4dca-9773-ecbc1e109c6a
md"""
## nanowire_tensorgrid
"""

# ╔═╡ ee5a4369-8ed2-44fd-b4f1-2a19f485abfb
grid_nanowire_tensorgrid=nanowire_tensorgrid(nrefs=2)

# ╔═╡ c3fe6feb-27b6-44e5-97a0-f0c3226f79cc
gridplot(  grid_nanowire_tensorgrid)

# ╔═╡ f78a3b26-a222-44f0-96b7-529e71dc0e81
md"""
## bimetal_tensorgrid
"""

# ╔═╡ 4cdc651c-54ea-4645-b8e1-0202c47e7f2c
grid_bimetal_tensorgrid=bimetal_tensorgrid()

# ╔═╡ cf0ca13b-b57b-4263-adf9-6d168ae75f0d
gridplot(grid_bimetal_tensorgrid)

# ╔═╡ baf477b5-a101-4a03-83a9-860540398a6e
md"""
## bimetal\_tensorgrid\_uniform
"""

# ╔═╡ 7dcb4899-5199-4201-82ee-64b046ea68ce
grid_bimetal_tensorgrid_uniform=bimetal_tensorgrid_uniform(nrefs=2)

# ╔═╡ Cell order:
# ╠═b8ec1d59-8a87-4cf0-bde8-790b420f71bc
# ╟─882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─b2d93b2f-78c9-4a4d-8cac-b04a625a3f75
# ╠═0b76e405-efd4-4c0b-8d53-486cf164e04c
# ╠═ca3fb071-ffda-46fc-ac22-0f3ecabca06a
# ╟─ef096445-8cf2-4d48-8011-e369527ad103
# ╠═d1c2cc0a-2d09-46b2-88b1-4c8358f42d2a
# ╟─e1d32a8e-49be-4c8e-9218-c01843368f73
# ╟─2608784c-d93e-41f6-a2b4-b2fcd63d8dd8
# ╠═2ede7b7b-1b03-49df-99fb-9137a877ab22
# ╠═5048badf-f655-4e8c-9455-c14a6913c53b
# ╟─3f648b46-dd2c-4ffd-9374-2b2493ac3d0c
# ╠═7c1178fc-b6db-48dd-83c0-cf8a266f2655
# ╠═7c54868c-756c-4a52-a07a-451d0677f6c1
# ╟─44a96db2-0faa-4536-a062-317d7e967f4f
# ╠═a786e949-aa71-475d-a82c-9190936b0f51
# ╠═719476ce-60c6-42ec-bf02-4679a0e67edc
# ╟─1c4624f1-d56c-47c5-869b-c29a4aab3f06
# ╠═ece590ef-61ad-47b6-9888-4e4c0d8163e7
# ╠═f85c2e83-d151-4647-a88e-15225c62631d
# ╟─0af4926e-eb19-4dca-9773-ecbc1e109c6a
# ╠═ee5a4369-8ed2-44fd-b4f1-2a19f485abfb
# ╠═c3fe6feb-27b6-44e5-97a0-f0c3226f79cc
# ╟─f78a3b26-a222-44f0-96b7-529e71dc0e81
# ╠═4cdc651c-54ea-4645-b8e1-0202c47e7f2c
# ╠═cf0ca13b-b57b-4263-adf9-6d168ae75f0d
# ╟─baf477b5-a101-4a03-83a9-860540398a6e
# ╠═7dcb4899-5199-4201-82ee-64b046ea68ce
