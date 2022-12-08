### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 6f8fc024-5e35-11ed-0828-8338123976ce
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
	using PyPlot
    default_plotter!(PlutoVista)
    TableOfContents()
end

# ╔═╡ 087da6a2-271f-45c7-98b5-a34a852ca8e7
md" ## Cross sections of different configurations"

# ╔═╡ 6cb742d3-519f-419f-8bc9-46b1c687a026
begin
		geometry = [30, 20, 10, 500]
	    scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
end

# ╔═╡ 5020c786-0560-4d00-a2ab-bca6b2c40fe7
let
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 1
	
	nrefs = 3
	refinement_width = 1/4*geometry[3]
	corner_refinement = false
	manual_refinement = true
	
    grid, cross_section = nanowire_tensorgrid_mirror(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ a0387376-89a9-4af1-945a-251ee6d668cd
let
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 2
	
	nrefs = 3
	refinement_width = 1/4*geometry[3]
	corner_refinement = false
	manual_refinement = true
	
    grid, cross_section = nanowire_tensorgrid_mirror(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ d60b1192-46be-4b6a-b48c-1144d6bf1257
let
	geometry = [15, 15, 10, 500]
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
	
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 3
	
	nrefs = 2
	refinement_width = 1/4*geometry[3]
	corner_refinement = false
	manual_refinement = true
	
    grid, cross_section = nanowire_tensorgrid_mirror(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ e3dd272a-59e9-4c3b-9e5c-ff0a061929e2
md" ## Further testing"

# ╔═╡ 4a50a0fa-0a2f-43ca-981f-b2239ec73885
md" #### Additional refinement at corners."

# ╔═╡ f9a24006-9153-4703-aa08-28277858746c
let
    grid, cross_section = nanowire_tensorgrid_mirror(; scale = scale,
        cut_levels = geometry[4]/2, nrefs = 1, refinement_width = nothing,
		shape = 2, corner_refinement = true, manual_refinement = false)
	@info num_nodes(grid)
	gridplot(cross_section)
end

# ╔═╡ 2ee98be7-3490-4f2a-b0e8-afa06494f3b4
md" #### Saving figure as a pdf."

# ╔═╡ cfd70770-9031-46cf-a0d6-321710b9899e
let
	grid, cross_section = nanowire_tensorgrid_mirror(; scale = scale,
        cut_levels = nothing, nrefs = 2, refinement_width = 1/2*geometry[3]/sqrt(3),
		shape = 2, corner_refinement = false, manual_refinement = true)
	vis = GridVisualizer(Plotter=PyPlot)
	gridplot!(vis,cross_section)
	PyPlot.savefig("cross_section.pdf")
	reveal(vis)
end

# ╔═╡ 1e3e6755-744a-4381-be41-05363664f335
md" #### Highlighting a point in cross section (useful in determining the different regions)."

# ╔═╡ 6e1c4ae9-15cc-436d-84b3-d9104bfd1b79
let
	shape = 2
	
	d1 = scale[1]
	d2 = scale[1] + scale[2]
	α = 1/2*geometry[3]/sqrt(3)
	δ = scale[3]

	grid, cross_section = nanowire_tensorgrid_mirror(; scale = scale,
		cut_levels = nothing, nrefs = 2, refinement_width = α, shape = shape,
        corner_refinement = false, manual_refinement = true)

    gridvis = GridVisualizer(Plotter=PyPlot,dim=2)
	gridplot!(gridvis,cross_section)
    ax = gridvis.subplots[1][:ax]
	# plot a thick red dot at given node
	if shape == 1
		ptx = -sqrt(3)/4*d2
		pty = -3/4*d2-(δ+α)/2 # -d2-δ+d2/4+(δ-α)/2
	elseif shape == 2
		ptx = -d1/2
		pty = -sqrt(3)/2*d2-(δ+α)/2
	end
    ax.scatter(ptx,pty,linewidths=0.1,color="red",zorder=2) 
    reveal(gridvis)
end

# ╔═╡ 65de36ff-5043-465f-beef-b2e2e7819ba2
begin
	md" ## Functions from NanoWiresJulia package
		(as per commit main/#a79de00af9758733804040aaf0cdef4715f6b8eb (08/12/22))
	"
end

# ╔═╡ 0658ec32-ba72-42ee-b8e5-6b2cd62d88b9
function asign_nodes!(p,builder,shape,d1,d2,δ)
	"""
		Assign nodes for different nanowire cross-section shapes.
	"""	
	if shape == 1
		p[1] = point!(builder,0,0)
		p[2] = point!(builder,0,d1)
		p[3] = point!(builder,-sqrt(3)/2*d1,d1/2)
		p[4] = point!(builder,-sqrt(3)/2*d1,-d1/2)
		p[5] = point!(builder,0,-d1)
		p[6] = point!(builder,0,d2)
		p[7] = point!(builder,-sqrt(3)/2*d2,d2/2)
		p[8] = point!(builder,-sqrt(3)/2*d2,-d2/2)
		p[9] = point!(builder,0,-d2)
		p[10] = point!(builder,0,-d2-δ)
		p[11] = point!(builder,-sqrt(3)/2*d2,-d2/2-δ)
	elseif shape == 2
	    p[1] = point!(builder,0,0)
	    p[2] = point!(builder,0,sqrt(3)/2*d1)
	    p[3] = point!(builder,-d1/2,sqrt(3)/2*d1)
	    p[4] = point!(builder,-d1,0)
	    p[5] = point!(builder,-d1/2,-sqrt(3)/2*d1)
	    p[6] = point!(builder,0,-sqrt(3)/2*d1)
	    p[7] = point!(builder,0,sqrt(3)/2*d2)
	    p[8] = point!(builder,-d2/2,sqrt(3)/2*d2)
	    p[9] = point!(builder,-d2,0)
	    p[10] = point!(builder,-d2/2,-sqrt(3)/2*d2)
	    p[11] = point!(builder,0,-sqrt(3)/2*d2)
		p[12] = point!(builder,0,-sqrt(3)/2*d2-δ)
	    p[13] = point!(builder,-d2/2,-sqrt(3)/2*d2-δ)
		p[14] = point!(builder,-d2,-δ)
	elseif shape == 3
	    p[1] = point!(builder,0,0)
	    p[2] = point!(builder,0,sqrt(3)/2*d1)
	    p[3] = point!(builder,-d1/2,sqrt(3)/2*d1)
	    p[4] = point!(builder,-d1,0)
	    p[5] = point!(builder,0,sqrt(3)/2*d2)
	    p[6] = point!(builder,-d2/2,sqrt(3)/2*d2)
	    p[7] = point!(builder,-d2,0)
		p[8] = point!(builder,0,sqrt(3)/2*(d2+δ))
		p[9] = point!(builder,-(d2+δ)/2,sqrt(3)/2*(d2+δ))
		p[10] = point!(builder,-d2-δ,0)
	end

	return p
end

# ╔═╡ 5966b4ec-d215-44b7-952c-1ea5643dd1ac
function interface_refinement!(p,builder,shape,d1,d2,δ,α)
	"""
		Allocate nodes at additonal regions along the material interface.
	"""
	if shape == 1
		if α > d2 - d1 || α > δ
			@warn "The refinement width cannot be larger than the shell's diameter or stressor width."
			return nothing
		end
		
		p[8] = point!(builder,-sqrt(3)/2*d2,-d2/2+α)
		p[9] = point!(builder,0,-d2+α)
		p[10] = point!(builder,0,-d2-δ)
		p[11] = point!(builder,-sqrt(3)/2*d2,-d2/2-δ)
		p[12] = point!(builder,-sqrt(3)/2*d2,-d2/2-α)
		p[13] = point!(builder,0,-d2-α)
		p[14] = point!(builder,0,-d2)
		p[15] = point!(builder,-sqrt(3)/2*d2,-d2/2)
	elseif shape == 2
		if α >= d2 - d1 || α >= δ/sqrt(3)
			@warn "The refinement width cannot be larger than the shell's diameter or δ/sqrt(3)."
			return nothing
		end
		
		p[9] = point!(builder,-d2+α/2,sqrt(3)/2*α)
		p[10] = point!(builder,(α-d2)/2,-sqrt(3)/2*(d2-α))
		p[11] = point!(builder,0,-sqrt(3)/2*(d2-α))
		p[15] = point!(builder,0,-sqrt(3)/2*d2)
		p[16] = point!(builder,-d2/2,-sqrt(3)/2*d2)
		p[17] = point!(builder,-d2,0)
		p[18] = point!(builder,-d2,-sqrt(3)*α)
		p[19] = point!(builder,-(d2+α)/2,-sqrt(3)/2*(d2+α))
		p[20] = point!(builder,0,-sqrt(3)/2*(d2+α))
	elseif shape == 3
		if α > d2 - d1 || α > δ
			@warn "The refinement width cannot be larger than the shell's or stressor's diameter."
			return nothing
		end
		p[5] = point!(builder,0,sqrt(3)/2*(d2-α))
		p[6] = point!(builder,(α-d2)/2,sqrt(3)/2*(d2-α))
		p[7] = point!(builder,-d2+α,0)
		p[8] = point!(builder,0,sqrt(3)/2*d2)
		p[9] = point!(builder,-d2/2,sqrt(3)/2*d2)
		p[10] = point!(builder,-d2,0)
		p[11] = point!(builder,0,sqrt(3)/2*(d2+α))
		p[12] = point!(builder,-(α+d2)/2,sqrt(3)/2*(d2+α))
		p[13] = point!(builder,-d2-α,0)
		p[14] = point!(builder,0,sqrt(3)/2*(d2+δ))
		p[15] = point!(builder,-(d2+δ)/2,sqrt(3)/2*(d2+δ))
		p[16] = point!(builder,-d2-δ,0)
	end
	
	return p
end

# ╔═╡ dd802a82-7e06-4ab0-9c60-bb7befec8500
function refine!(builder,shape,d1,d2,δ,α)
	"""
		Assign points along the middle of the each interface region to enable a
		uniform refinement.
	"""
	if shape == 1
		num_pts = trunc(Int, 4*δ/α)
		for n = 0 : num_pts
			g = n/num_pts
			# convex combination between p8 & p15 midpoint and p9 & p14 midpoint
			px = g*(-sqrt(3)/2*d2)
			py = g*(-d2/2+α/2) + (1-g)*(-d2+α/2)
			point!(builder,px,py)
			# convex combination between p15 & p12 midpoint and p14 & p13 midpoint
			px = g*(-sqrt(3)/2*d2)
			py = g*(-d2/2-α/2) + (1-g)*(-d2-α/2)
			point!(builder,px,py)
		end
	elseif shape == 2
		num_pts = trunc(Int, 4*δ/(sqrt(3)*α))
		for n = 0 : num_pts
			g = n/num_pts
			## adding points along the horizontal interface
			# convex combination between p10 & p16 midpoint and p11 & p15 midpoint
			px = g*(α/2-d2)/2
			py = -sqrt(3)/2*(d2-α/2)
			point!(builder,px,py)
			# convex combination between p16 & p19 midpoint and p15 & p20 midpoint
			px = g*(-(d2+α/2)/2)
			py = -sqrt(3)/2*(d2+α/2)
			point!(builder,px,py)
		end
		num_pts = 2*num_pts
		for n = 0 : num_pts-1
			g = n/num_pts
			## adding points along the left interface
			# convex combination between p9 & p17 midpoint and p10 & p16 midpoint
			px = g*(-d2+α/4) + (1-g)*(α/2-d2)/2
			py = g*sqrt(3)/2*α/2 + (1-g)*(-sqrt(3)/2*(d2-α/2))
			point!(builder,px,py)
			# convex combination between p17 & p18 midpoint and p16 & p19 midpoint
			px = g*(-d2) + (1-g)*(-(d2+α/2)/2)
			py = g*(-sqrt(3)*α/2) + (1-g)*(-sqrt(3)/2*(d2+α/2))
			point!(builder,px,py)
		end
	elseif shape == 3
		num_pts = trunc(Int, 2*δ/α)
		for n = 0 : num_pts-1
			g = n/num_pts
			# convex combination between p5 & p8 midpoint and p6 & p9 midpoint
			px = (1-g)*(α/2-d2)/2
			py = sqrt(3)/2*(d2-α/2)
			point!(builder,px,py)
			# convex combination between p8 & p11 midpoint and p9 & p12 midpoint
			px = (1-g)*(-α/2-d2)/2
			py = sqrt(3)/2*(d2+α/2)
			point!(builder,px,py)
		end
		num_pts = 2*num_pts
		for n = 0 : num_pts-1
			g = n/num_pts
			# convex combination between p6 & p9 midpoint and p7 & p10 midpoint
			px = g*(α/2-d2)/2 + (1-g)*(-d2+α/2)
			py = g*sqrt(3)/2*(d2-α/2)
			point!(builder,px,py)
			# convex combination between p9 & p12 midpoint and p10 & p13 midpoint
			px = g*(-α/2-d2)/2 + (1-g)*(-d2-α/2)
			py = g*sqrt(3)/2*(d2+α/2)
			point!(builder,px,py)
		end
	end

end

# ╔═╡ e502a079-e1f9-4484-9124-ecc9cf49363c
function assign_core_shell_edges!(p,builder,shape)
	"""
		Assign edges.
	"""
	if shape == 1
		facetregion!(builder,1) # core region
		facet!(builder,p[1],p[2])
		facet!(builder,p[2],p[3])
		facet!(builder,p[3],p[4])
		facet!(builder,p[4],p[5])
		facet!(builder,p[5],p[1])

		facetregion!(builder,2) # shell region
		facet!(builder,p[2],p[6])
		facet!(builder,p[6],p[7])
		facet!(builder,p[7],p[8])
		facet!(builder,p[8],p[9])
		facet!(builder,p[9],p[5])
	elseif shape == 2
		facetregion!(builder,1) # core region
		facet!(builder,p[1],p[2])
		facet!(builder,p[2],p[3])
		facet!(builder,p[3],p[4])
		facet!(builder,p[4],p[5])
		facet!(builder,p[5],p[6])
		facet!(builder,p[6],p[1])
	
		facetregion!(builder,2) # shell region
		facet!(builder,p[2],p[7])
		facet!(builder,p[7],p[8])
		facet!(builder,p[8],p[9])
		facet!(builder,p[9],p[10])
		facet!(builder,p[10],p[11])
		facet!(builder,p[11],p[6])
	elseif shape == 3
		facetregion!(builder,1) # core region
		facet!(builder,p[1],p[2])
		facet!(builder,p[2],p[3])
		facet!(builder,p[3],p[4])
		facet!(builder,p[4],p[1])
	
		facetregion!(builder,2) # shell region
		facet!(builder,p[2],p[5])
		facet!(builder,p[5],p[6])
		facet!(builder,p[6],p[7])
		facet!(builder,p[7],p[4])
	end
end

# ╔═╡ db067db2-2b9a-4bd4-ac8a-abd58fb24094
function assign_stressor_edges!(p,builder,shape,refinement_width)
	"""
		Assign stressor edges.
	"""
	if refinement_width == nothing
		facetregion!(builder,3) # stressor region
		if shape == 1
			facet!(builder,p[9],p[10])
			facet!(builder,p[10],p[11])
			facet!(builder,p[11],p[8])
		elseif shape == 2
			facet!(builder,p[11],p[12])
			facet!(builder,p[12],p[13])
			facet!(builder,p[13],p[14])
			facet!(builder,p[14],p[9])
		elseif shape == 3
			facet!(builder,p[5],p[8])
			facet!(builder,p[8],p[9])
			facet!(builder,p[9],p[10])
			facet!(builder,p[10],p[7])
		end
	else
		if shape == 1
			facetregion!(builder,2) # interface shell region
			facet!(builder,p[9],p[14])
			facet!(builder,p[14],p[15])
			facet!(builder,p[15],p[8])

			facetregion!(builder,3) # interface stressor region
			facet!(builder,p[15],p[12])
			facet!(builder,p[12],p[13])
			facet!(builder,p[13],p[14])
		
			facetregion!(builder,3) # stressor region
			facet!(builder,p[13],p[10])
			facet!(builder,p[10],p[11])
			facet!(builder,p[11],p[12])
		elseif shape == 2
			facetregion!(builder,2) # interface shell region
			facet!(builder,p[11],p[15])
			facet!(builder,p[15],p[16])
			facet!(builder,p[16],p[17])
			facet!(builder,p[17],p[9])
		
			facetregion!(builder,3) # interface stressor region
			facet!(builder,p[17],p[18])
			facet!(builder,p[18],p[19])
			facet!(builder,p[19],p[20])
			facet!(builder,p[20],p[15])
		
			facetregion!(builder,3) # stressor region
			facet!(builder,p[20],p[12])
			facet!(builder,p[12],p[13])
			facet!(builder,p[13],p[14])
			facet!(builder,p[14],p[18])
		elseif shape == 3
			facetregion!(builder,2) # interface shell region
			facet!(builder,p[5],p[8])
			facet!(builder,p[8],p[9])
			facet!(builder,p[9],p[10])
			facet!(builder,p[10],p[7])
		
			facetregion!(builder,3) # interface stressor region
			facet!(builder,p[8],p[11])
			facet!(builder,p[11],p[12])
			facet!(builder,p[12],p[13])
			facet!(builder,p[13],p[10])
		
			facetregion!(builder,3) # stressor region
			facet!(builder,p[11],p[14])
			facet!(builder,p[14],p[15])
			facet!(builder,p[15],p[16])
			facet!(builder,p[16],p[13])
		end
	end
end

# ╔═╡ 50565fe6-d590-42bd-90d9-5297c76df062
function assign_regions!(p,builder,shape,refinement_width,nrefs,d1,d2,δ)
	"""
		Assign regions.
	"""
	A_core = 1/2 * (3*sqrt(3)/2*d1^2)
    A_shell = 1/2 * (3*sqrt(3)/2*(d2^2 - d1^2))

	vol_factor_core = 4.0^-nrefs
	vol_factor_shell = 4.0^-nrefs
	vol_factor_stressor = 4.0^-nrefs

	if shape == 1

		A_stressor = δ*d2

		if refinement_width == nothing
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-sqrt(3)/4*d2,-d2))
		else
			α = refinement_width
			A_interface_interior = d2*α
			# the refinement area of the stressor region is the same as the
			# refinement region in the shell
			A_interface_exterior = A_interface_interior
			A_shell = A_shell - A_interface_interior
			A_stressor = A_stressor - A_interface_exterior
	
			vol_factor_core = 4.0^-(nrefs-1)
			vol_factor_interface = 4.0^-nrefs
			vol_factor_stressor = 4.0^-(nrefs-1)
		
			cellregion!(builder,2) # material 2 (interface shell)
			maxvolume!(builder,A_interface_interior*vol_factor_interface)
			regionpoint!(builder,(-sqrt(3)/4*d2,1/2*(-3/2*d2+α)))
		
			cellregion!(builder,3) # material 3 (interface stressor)
			maxvolume!(builder,A_interface_exterior*vol_factor_interface)
			regionpoint!(builder,(-sqrt(3)/4*d2,1/2*(-3/2*d2-α)))
			
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-sqrt(3)/4*d2,-3/4*d2-(δ+α)/2)) # -d2-δ+d2/4+(δ-α)/2
		end
	
		cellregion!(builder,1) # material 1 (core)
		maxvolume!(builder,A_core*vol_factor_core)
		regionpoint!(builder,(-sqrt(3)/4*d1,0))
	
		cellregion!(builder,2) # material 2 (shell)
		maxvolume!(builder,A_shell*vol_factor_shell)
		regionpoint!(builder,(-sqrt(3)/2*(d1+d2)/2),0)

	elseif shape == 2

		A_stressor = 1/2 * (δ*(d2 + δ/sqrt(3)) + δ*d2)

		if refinement_width == nothing
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-d1/2,-sqrt(3)/2*d2-δ/2))
		else
			α = refinement_width
			A_interface_interior = (2*d2-α)*sqrt(3)/2*α + d2*α
			# the refinement area of the stressor region is almost the same as the
			# refinement region in the shell
			A_interface_exterior = A_interface_interior
			A_shell = A_shell - A_interface_interior
			A_stressor = A_stressor - A_interface_exterior
	
			vol_factor_core = 4.0^-(nrefs-1)
			vol_factor_interface = 4.0^-nrefs
			vol_factor_stressor = 4.0^-(nrefs-1)

			cellregion!(builder,2) # material 2 (interface shell)
			maxvolume!(builder,A_interface_interior*vol_factor_interface)
			regionpoint!(builder,((α/2-d2)/2,-sqrt(3)/2*(d2-α/2)))
		
			cellregion!(builder,3) # material 3 (interface stressor)
			maxvolume!(builder,A_interface_exterior*vol_factor_interface)
			regionpoint!(builder,(-(d2+α/2)/2,-sqrt(3)/2*(d2+α/2)))
			
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-d1/2,-sqrt(3)/2*d2-(δ+α)/2))
		end
	
		cellregion!(builder,1) # material 1 (core)
		maxvolume!(builder,A_core*vol_factor_core)
		regionpoint!(builder,(-d1/2,0))
	
		cellregion!(builder,2) # material 2 (shell)
		maxvolume!(builder,A_shell*vol_factor_shell)
		regionpoint!(builder,(-d1/2,sqrt(3)/2*(d1+d2)/2))

	elseif shape == 3

		A_core /= 2
		A_shell /= 2
		A_stressor = 1/4 * (3*sqrt(3)/2*((d2+δ)^2 - d2^2))

		if refinement_width == nothing
			cellregion!(builder,2) # material 2 (shell)
			maxvolume!(builder,A_shell*vol_factor_shell)
			regionpoint!(builder,(-d2/4,sqrt(3)/2*(d1+d2)/2))
			
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-(d2+δ)/4,sqrt(3)/2*(2*d2+δ)/2))
		else
			α = refinement_width
			A_interface_interior = 3/2 * (2*d2-α)*sqrt(3)/2*α
			A_interface_exterior = 3/2 * (2*d2+α)*sqrt(3)/2*α
			A_shell = A_shell - A_interface_interior
			A_stressor = A_stressor - A_interface_exterior
	
			vol_factor_core = 4.0^-(nrefs-1)
			vol_factor_shell = 4.0^-(nrefs-1)
			vol_factor_interface = 4.0^-nrefs
			vol_factor_stressor = 4.0^-(nrefs-1)

			cellregion!(builder,2) # material 2 (shell)
			maxvolume!(builder,A_shell*vol_factor_shell)
			regionpoint!(builder,(-d2/4,sqrt(3)/2*(d1+d2-α)/2))
			
			cellregion!(builder,2) # material 2 (interface shell)
			maxvolume!(builder,A_interface_interior*vol_factor_interface)
			regionpoint!(builder,(-d2/4,sqrt(3)/2*(d2-α/2)))
		
			cellregion!(builder,3) # material 3 (interface stressor)
			maxvolume!(builder,A_interface_exterior*vol_factor_interface)
			regionpoint!(builder,(-d2/4,sqrt(3)/2*(d2+α/2)))
			
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-(d2+δ)/4,sqrt(3)/2*(2*d2+δ+α)/2))
		end
		
		cellregion!(builder,1) # material 1 (core)
		maxvolume!(builder,A_core*vol_factor_core)
		regionpoint!(builder,(-d1/4,sqrt(3)/2*d1/2))

	end

end

# ╔═╡ 5bb8e9ca-8669-46a1-97d4-72a912465e9e
function nanowire_tensorgrid_mirror!(; scale = [1,1,1,1], shape = 1,
	nrefs = 1, z_nrefs = 2, z_levels_dist = 100, cut_levels = scale[4]/2,
	refinement_width = nothing, corner_refinement = false, manual_refinement = false, max_nodes = 20)
    
    @info "Generating nanowire grid for scale = $scale"

    builder = SimplexGridBuilder(Generator=Triangulate)
	p::Array{Int64,1} = zeros(Int64,max_nodes)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ  = scale[3]
	α = refinement_width    

    ## assign nodes
	p = asign_nodes(p,builder,shape,d1,d2,δ)

	## assign extra nodes for refinement accross the material interface 
	if refinement_width !== nothing
		p = interface_refinement(p,builder,shape,d1,d2,δ,refinement_width)
		
		if manual_refinement == true
			refine(builder,shape,d1,d2,δ,refinement_width)
		end
	end

	## building edges
	assign_core_shell_edges(p,builder,shape)
	assign_stressor_edges(p,builder,shape,refinement_width)

	## assigning regions
	assign_regions(p,builder,shape,refinement_width,nrefs,d1,d2,δ)
	
	# optional refinement at interface corners 
	if corner_refinement == true
        function unsuitable(x1,y1,x2,y2,x3,y3, area)
            bary = [(x1+x2+x3)/3,(y2+y2+y3)/3]
			dist = min(norm(bary-refinement_center1),norm(bary-refinement_center2))
			if shape == 3
				dist = min(dist,norm(bary-refinement_center3))
			end
            if area > 1.5*dist
                return 1
            else
                return 0
            end
        end

		if shape == 1
			refinement_center1 = [-sqrt(3)/2*d2,-d2/2]
        	refinement_center2 = [0,-d2]
		else
	        refinement_center1 = [-d2,0]
    	    refinement_center2 = [-d2/2,-sqrt(3)/2*d2]
			refinement_center3 = [-d2/2,sqrt(3)/2*d2]
		end
		options!(builder, unsuitable=unsuitable)
	end

	# build grid
	xgrid = simplexgrid(builder)
	
	# mirror grid
	xgrid_flipped=deepcopy(xgrid)
    xgrid_flipped[Coordinates][1,:] .*= -1
    xgrid=glue(xgrid,xgrid_flipped; interface=4)
    bfacemask!(xgrid,[0,-(d2+δ)],[0,d2+δ],0)
	if shape == 3
		xgrid_flipped=deepcopy(xgrid)
    	xgrid_flipped[Coordinates][2,:] .*= -1
    	xgrid=glue(xgrid,xgrid_flipped; interface=4)
    	bfacemask!(xgrid,[-(d2+δ),0],[d2+δ,0],0)
	end
	xgrid_cross_section=deepcopy(xgrid)

	# tensor grid in the z-direction
	hz_factor = 2.0^-z_nrefs
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
	
	# generating final grid 
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

	return xgrid, xgrid_cross_section
end

# ╔═╡ Cell order:
# ╠═6f8fc024-5e35-11ed-0828-8338123976ce
# ╟─087da6a2-271f-45c7-98b5-a34a852ca8e7
# ╠═6cb742d3-519f-419f-8bc9-46b1c687a026
# ╠═5020c786-0560-4d00-a2ab-bca6b2c40fe7
# ╠═a0387376-89a9-4af1-945a-251ee6d668cd
# ╠═d60b1192-46be-4b6a-b48c-1144d6bf1257
# ╟─e3dd272a-59e9-4c3b-9e5c-ff0a061929e2
# ╟─4a50a0fa-0a2f-43ca-981f-b2239ec73885
# ╠═f9a24006-9153-4703-aa08-28277858746c
# ╟─2ee98be7-3490-4f2a-b0e8-afa06494f3b4
# ╠═cfd70770-9031-46cf-a0d6-321710b9899e
# ╟─1e3e6755-744a-4381-be41-05363664f335
# ╠═6e1c4ae9-15cc-436d-84b3-d9104bfd1b79
# ╟─65de36ff-5043-465f-beef-b2e2e7819ba2
# ╟─0658ec32-ba72-42ee-b8e5-6b2cd62d88b9
# ╟─5966b4ec-d215-44b7-952c-1ea5643dd1ac
# ╟─dd802a82-7e06-4ab0-9c60-bb7befec8500
# ╟─e502a079-e1f9-4484-9124-ecc9cf49363c
# ╟─db067db2-2b9a-4bd4-ac8a-abd58fb24094
# ╟─50565fe6-d590-42bd-90d9-5297c76df062
# ╟─5bb8e9ca-8669-46a1-97d4-72a912465e9e
