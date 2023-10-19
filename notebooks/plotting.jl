### A Pluto.jl notebook ###
# v0.19.29

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

# ╔═╡ 902d154e-5a51-4955-b3c7-9eda226a7125
md" # Nanowire"

# ╔═╡ 087da6a2-271f-45c7-98b5-a34a852ca8e7
md" ## Cross sections of different configurations"

# ╔═╡ 6cb742d3-519f-419f-8bc9-46b1c687a026
begin
		geometry = [30, 20, 10, 500]
	    scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
end

# ╔═╡ c9bedd54-b2f0-416e-ae9e-f3e83e9f7d76
let
	dbulk = 50
	dstressor = 10
	dcore = 30
	geometry = [dcore, dbulk-dcore, dstressor, 500]
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 2
	
	nrefs = 0
	refinement_width = nothing
	corner_refinement = false
	manual_refinement = false

	grid = nanowire_grid(; scale=scale, reflevel=nrefs)
	gridplot!(vis[1,1],grid)
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

# ╔═╡ 1e3e6755-744a-4381-be41-05363664f335
md" #### Highlighting a point in cross section. Useful to determine regions."

# ╔═╡ 65de36ff-5043-465f-beef-b2e2e7819ba2
md" ## Functions from NanoWiresJulia package"

# ╔═╡ 0658ec32-ba72-42ee-b8e5-6b2cd62d88b9
function asign_nodes(p,builder,shape,d1,d2,δ)
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
		p[8] = point!(builder,0,sqrt(3)/2*d2+δ)
		p[9] = point!(builder,-(d2+2/sqrt(3)*δ)/2,sqrt(3)/2*d2+δ)
		p[10] = point!(builder,-d2-2/sqrt(3)*δ,0)
	end

	return p
end

# ╔═╡ 5966b4ec-d215-44b7-952c-1ea5643dd1ac
function interface_refinement(p,builder,shape,d1,d2,δ,α)
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
		p[14] = point!(builder,0,sqrt(3)/2*d2+δ)
		p[15] = point!(builder,-(d2+2/sqrt(3)*δ)/2,sqrt(3)/2*d2+δ)
		p[16] = point!(builder,-d2-2/sqrt(3)*δ,0)
	end
	
	return p
end

# ╔═╡ dd802a82-7e06-4ab0-9c60-bb7befec8500
function refine(builder,shape,d1,d2,δ,α)
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
		num_pts = 10 #min(trunc(Int, 4*δ/(sqrt(3)*α)),10)
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
		num_pts = trunc(Int, d2/4)+2
		for n = 0 : num_pts
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
function assign_core_shell_edges(p,builder,shape)
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
function assign_stressor_edges(p,builder,shape,refinement_width)
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
function assign_regions(p,builder,shape,refinement_width,nrefs,d1,d2,δ)
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
		A_stressor = 1/4 * (3*sqrt(3)/2*((d2+2/sqrt(3)*δ)^2 - d2^2))

		if refinement_width == nothing
			cellregion!(builder,2) # material 2 (shell)
			maxvolume!(builder,A_shell*vol_factor_shell)
			regionpoint!(builder,(-d2/4,sqrt(3)/2*(d1+d2)/2))
			
			cellregion!(builder,3) # material 3 (stressor)
			maxvolume!(builder,A_stressor*vol_factor_stressor)
			regionpoint!(builder,(-(d2+2/sqrt(3)*δ)/4,sqrt(3)/2*d2+δ/2))
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
			regionpoint!(builder,(-(d2+2/sqrt(3)*δ)/4,sqrt(3)/2*(2*d2+α)/2+δ/2))
		end
		
		cellregion!(builder,1) # material 1 (core)
		maxvolume!(builder,A_core*vol_factor_core)
		regionpoint!(builder,(-d1/4,sqrt(3)/2*d1/2))

	end

end

# ╔═╡ 5bb8e9ca-8669-46a1-97d4-72a912465e9e
function nanowire_tensorgrid_mirror!(; scale = [1,1,1,1], shape = 1,
	nrefs = 1, z_nrefs = 2, z_levels_dist = 100, cut_levels = scale[4]/2,
	refinement_width = nothing, corner_refinement = false, manual_refinement = false,
	rotate = true, max_nodes = 20)
    
    @info "Generating nanowire grid for scale = $scale"

    builder = SimplexGridBuilder(Generator=Triangulate)
	p::Array{Int64,1} = zeros(Int64,max_nodes)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ  = scale[3]
	
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
	xgrid_flipped = deepcopy(xgrid)
	xgrid_flipped[Coordinates][1,:] .*= -1
	xgrid = glue(xgrid,xgrid_flipped; interface=4)
    bfacemask!(xgrid,[0,-(d2+δ)],[0,d2+δ],0)
	if shape == 3
		xgrid_flipped = deepcopy(xgrid)
    	xgrid_flipped[Coordinates][2,:] .*= -1
    	xgrid = glue(xgrid,xgrid_flipped; interface=4)
    	bfacemask!(xgrid,[-(d2+2/sqrt(3)*δ),0],[d2+2/sqrt(3)*δ,0],0)
	end
	
	# clockwise grid rotation by angle θ
	θ = rotate * pi/180
	xgrid_rotated = deepcopy(xgrid)
	xgrid_rotated[Coordinates][1,:] = cos(θ)*xgrid[Coordinates][1,:] + 
		sin(θ)*xgrid[Coordinates][2,:]
	xgrid_rotated[Coordinates][2,:] = -sin(θ)*xgrid[Coordinates][1,:] + 
		cos(θ)*xgrid[Coordinates][2,:]
	xgrid = deepcopy(xgrid_rotated)
	repair_grid!(xgrid)
	
	# save cross-section grid
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

# ╔═╡ 6a1e0b51-cf24-4f24-aaf0-39fe3beb8dcd
let
	dbulk = 50
	dstressor = 10
	dcore = 30
	geometry = [dcore, dbulk-dcore, dstressor, 2000]
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 2
	
	nrefs = 0
	refinement_width = nothing
	corner_refinement = false
	manual_refinement = false
	rotate = 0

    grid, cross_section = nanowire_tensorgrid_mirror!(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
		corner_refinement=corner_refinement, manual_refinement=manual_refinement,
		rotate=rotate)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ 5020c786-0560-4d00-a2ab-bca6b2c40fe7
let
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 1
	
	nrefs = 3
	refinement_width = 1/4*geometry[3]
	corner_refinement = false
	manual_refinement = true
	rotate = 45
	
    grid, cross_section = nanowire_tensorgrid_mirror!(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement,
		rotate=rotate)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ a0387376-89a9-4af1-945a-251ee6d668cd
let
	dbulk = 50
	dstressor = 10
	dcore = 20
	geometry = [dcore, dbulk-dcore, dstressor, 500]
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 2
	
	nrefs = 3
	refinement_width = (dstressor >= 10 ? 4 : dbulk/25+1)/sqrt(3)
	corner_refinement = false
	manual_refinement = true
	rotate = 45
	
    grid, cross_section = nanowire_tensorgrid_mirror!(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement,
		rotate=rotate)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ f478eda3-6192-4f69-ba39-08bebeaa7847
let
	geometry = [30, 20, 15, 500]
	
	scale = [geometry[1],geometry[2],geometry[3]*2/sqrt(3),geometry[4]]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 3
	
	nrefs = 2
	refinement_width = 1/4*geometry[2]
	corner_refinement = false
	manual_refinement = true
	rotate = 0
	
    grid, cross_section = nanowire_tensorgrid_mirror!(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement,
		rotate=rotate)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ d60b1192-46be-4b6a-b48c-1144d6bf1257
let
	geometry = [7.5, 7.5, 50, 500]
	
	scale = [geometry[1],geometry[2],geometry[3]*sqrt(3)/2,geometry[4]]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 3
	
	nrefs = 2
	refinement_width = nothing #1/2*geometry[2]
	corner_refinement = false
	manual_refinement = false
	rotate = 0
	
    grid, cross_section = nanowire_tensorgrid_mirror!(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement,
		rotate = false)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ ecf47d9b-efbb-467a-b924-3195a6cdd9cf
let
	geometry = [15*sqrt(3), 15*sqrt(3), 50, 500]
	
	scale = [geometry[1]/sqrt(3),geometry[2]/sqrt(3),geometry[3],geometry[4]]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	shape = 3
	
	nrefs = 0
	refinement_width = nothing #1/4*geometry[2]
	corner_refinement = false
	manual_refinement = false
	rotate = 0
	
    grid, cross_section = nanowire_tensorgrid_mirror!(; scale=scale, shape=shape,
        cut_levels=nothing, nrefs=nrefs, refinement_width=refinement_width,
        corner_refinement=corner_refinement, manual_refinement=manual_refinement,
		rotate=rotate)
	gridplot!(vis[1,1],cross_section)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ cfd70770-9031-46cf-a0d6-321710b9899e
let
	grid, cross_section = nanowire_tensorgrid_mirror!(; scale = scale,
        cut_levels = nothing, nrefs = 2, refinement_width = 1/2*geometry[3]/sqrt(3),
		shape = 2, corner_refinement = false, manual_refinement = true, rotate = 0)
	vis = GridVisualizer(Plotter=PyPlot)
	gridplot!(vis,cross_section)
	PyPlot.savefig("cross_section.pdf")
	reveal(vis)
end

# ╔═╡ 6e1c4ae9-15cc-436d-84b3-d9104bfd1b79
let
	shape = 3
	rotate = 90
	
	d1 = scale[1]
	d2 = scale[1] + scale[2]
	α = 1/2*geometry[3]/sqrt(3)
	δ = scale[3]

	grid, cross_section = nanowire_tensorgrid_mirror!(; scale = scale,
		cut_levels = nothing, nrefs = 2, refinement_width = α, shape = shape,
        corner_refinement = false, manual_refinement = true, rotate = rotate)

    gridvis = GridVisualizer(Plotter=PyPlot,dim=2)
	gridplot!(gridvis,cross_section)
    ax = gridvis.subplots[1][:ax]
	# plot a thick red dot at given node
	if shape == 1
		ptx = -3/4*d2-(δ+α)/2 # -d2-δ+d2/4+(δ-α)/2
		pty = sqrt(3)/4*d2
		ptx = -(geometry[1]+geometry[2])/sqrt(3)
		pty = 0
		@info ptx,pty
	elseif shape == 2 || shape == 3 
		ptx = -sqrt(3)/2*d2-(δ+α)/2
		pty = d1/2
		ptx = -(geometry[1]+geometry[2])/2
		pty = 0
		@info ptx,pty
	end
    ax.scatter(ptx,pty,linewidths=0.1,color="red",zorder=2) 
    reveal(gridvis)
end

# ╔═╡ 4a50c199-2d06-420a-b797-3faa40ab14b8
md" # Bimetal"

# ╔═╡ 1e6e16d0-a4e4-43e0-a67b-88656c1f5247
let
	scale = [50,100,500]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	nrefs = 1
	mb = 0.5
	hz = 50
	
    grid, xgrid_cross_section = bimetal_tensorgrid_uniform(; scale=scale, nrefs=nrefs, material_border=mb, hz=hz)
	gridplot!(vis[1,1],grid)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ f3203adb-6624-4cb8-bab7-825e9e2d956a
function bimetal_tensorgrid_uniform1(; scale = [1,1,1], nrefs = 1, material_border = 0.5)

    @info "Generating bimetal 3D grid for scale = $scale and middle interface at $material_border of height $(scale[2])"

    W = scale[1]
    H = scale[2]
    Z = scale[3]
    h1 = round(scale[2]*material_border)
    h2 = round(scale[2]*(1 - material_border))
    hz_factor = 2.0^-nrefs

    hx = W/4*2.0^-nrefs
    hy = min(h1,h2)*2.0^-nrefs
    hz = 50 * hz_factor

    XX = 0:hx:W
    #YY = Array{Float64,1}(0:(h1-hy)/2:h1-hy)
    #append!(YY, Array{Float64,1}(LinRange(h1-hy,h1+hy,3)[2:end-1]))
    #append!(YY, Array{Float64,1}(h1+hy:(h2-hy)/2:H))
    YY = 0:hy:H
    ZZ = Array{Float64,1}(0:hz:Z)

    xgrid = simplexgrid(XX,YY)
	
	# assigning region numbers: core region = 1, stressor region = 2
    cellmask!(xgrid,[0,0],[W,h1],1)
    cellmask!(xgrid,[0,h1],[W,H],2)
	
	xgrid_cross_section = deepcopy(xgrid)
    xgrid = simplexgrid(xgrid, ZZ)
    xgrid = uniform_refine(xgrid,nrefs)
    # the offsets lead to the following boundary regions:
    # 1 - 6 = sides core & bottom
    # 7 = bottom core
    # 8 = bottom stressor
    # 9 = top core
    # 10 = top stressor

    # boundary faces
    bfacemask!(xgrid,[0,0,0],[W,h1,Z],1)  # side core left
    bfacemask!(xgrid,[0,h1,0],[0,H,Z],2)  # side stressor left
    bfacemask!(xgrid,[0,H,0],[W,H,Z],3)   # side stressor
    bfacemask!(xgrid,[W,0,0],[W,h1,Z],4)  # side core right
    bfacemask!(xgrid,[W,h1,0],[W,H,Z],5)  # side stressor right
    bfacemask!(xgrid,[0,0,0],[W,0,Z],6)   # side core
    bfacemask!(xgrid,[0,0,0],[W,h1,0],7)  # bottom core
    bfacemask!(xgrid,[0,h1,0],[W,H,0],8)  # bottom stressor
    bfacemask!(xgrid,[0,0,Z],[W,h1,Z],9)  # top core
    bfacemask!(xgrid,[0,h1,Z],[W,H,Z],10) # top stressor

    return xgrid, xgrid_cross_section

end

# ╔═╡ 8f3027d3-fb5f-4f24-ab6b-4dc4526a4f99
let
	scale = [200,50,2000]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	nrefs = 1
	mb = 0.5
	
    grid,xgrid_cross_section = bimetal_tensorgrid_uniform1(; scale=scale, nrefs=nrefs, material_border=mb)
	gridplot!(vis[1,1],grid)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ 79ded2f2-eb09-4114-abda-d79b6111fbb5
function bimetal_tensorgrid_uniform2(; scale = [1,1,1], nrefs = 1, material_border = 0.5)

    @info "Generating bimetal 3D grid for scale = $scale and middle interface at $material_border of height $(scale[2])"

    W = scale[1]
    H = scale[2]
    Z = scale[3]
    h1 = round(scale[2]*material_border)
    h2 = round(scale[2]*(1 - material_border))
    hz_factor = 2.0^-nrefs

    hx = W/2 # W/4*2.0^-nrefs
    hy = H/2 # min(h1,h2)/2*2.0^-nrefs
    hz = 50 # 100 * hz_factor

    XX = 0:hx:W
    #YY = Array{Float64,1}(0:(h1-hy)/2:h1-hy)
    #append!(YY, Array{Float64,1}(LinRange(h1-hy,h1+hy,3)[2:end-1]))
    #append!(YY, Array{Float64,1}(h1+hy:(h2-hy)/2:H))
    YY = 0:hy:H
    ZZ = Array{Float64,1}(0:hz:Z)

    xgrid = simplexgrid(YY,XX)

	# assigning region numbers: core region = 1, stressor region = 2
    cellmask!(xgrid,[0,0],[h1,W],1)
    cellmask!(xgrid,[h1,0],[H,W],2)
	
    xgrid_cross_section = deepcopy(xgrid)
	coords = deepcopy(xgrid[Coordinates])
	xgrid[Coordinates][1,:] = -coords[2,:]
	xgrid[Coordinates][2,:] = coords[1,:]
	xgrid[Coordinates][1,:] .+= W
    xgrid = simplexgrid(xgrid, ZZ)
    xgrid = uniform_refine(xgrid,nrefs)
    # the offsets lead to the following boundary regions:
    # 1 - 6 = sides core & bottom
    # 7 = bottom core
    # 8 = bottom stressor
    # 9 = top core
    # 10 = top stressor

    # boundary faces
    bfacemask!(xgrid,[0,0,0],[W,h1,Z],1)  # side core left
    bfacemask!(xgrid,[0,h1,0],[0,H,Z],2)  # side stressor left
    bfacemask!(xgrid,[0,H,0],[W,H,Z],3)   # side stressor
    bfacemask!(xgrid,[W,0,0],[W,h1,Z],4)  # side core right
    bfacemask!(xgrid,[W,h1,0],[W,H,Z],5)  # side stressor right
    bfacemask!(xgrid,[0,0,0],[W,0,Z],6)   # side core
    bfacemask!(xgrid,[0,0,0],[W,h1,0],7)  # bottom core
    bfacemask!(xgrid,[0,h1,0],[W,H,0],8)  # bottom stressor
    bfacemask!(xgrid,[0,0,Z],[W,h1,Z],9)  # top core
    bfacemask!(xgrid,[0,h1,Z],[W,H,Z],10) # top stressor

    return xgrid, xgrid_cross_section

end

# ╔═╡ 4664553e-d37a-41f7-a439-c833e8ba95a9
let
	scale = [20,100,200]
	vis = GridVisualizer(Plotter=PlutoVista,layout=(1,1),resolution=(700,700))

	nrefs = 1
	mb = 0.5
	hz = 50
	
    grid, xgrid_cross_section = bimetal_tensorgrid_uniform(; scale=scale, nrefs=nrefs, material_border=mb, hz=hz)
	gridplot!(vis[1,1],grid)
	@info num_nodes(grid)
	reveal(vis)
end

# ╔═╡ Cell order:
# ╠═6f8fc024-5e35-11ed-0828-8338123976ce
# ╟─902d154e-5a51-4955-b3c7-9eda226a7125
# ╟─087da6a2-271f-45c7-98b5-a34a852ca8e7
# ╠═6cb742d3-519f-419f-8bc9-46b1c687a026
# ╠═c9bedd54-b2f0-416e-ae9e-f3e83e9f7d76
# ╠═6a1e0b51-cf24-4f24-aaf0-39fe3beb8dcd
# ╠═5020c786-0560-4d00-a2ab-bca6b2c40fe7
# ╠═a0387376-89a9-4af1-945a-251ee6d668cd
# ╠═f478eda3-6192-4f69-ba39-08bebeaa7847
# ╠═d60b1192-46be-4b6a-b48c-1144d6bf1257
# ╠═ecf47d9b-efbb-467a-b924-3195a6cdd9cf
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
# ╠═5bb8e9ca-8669-46a1-97d4-72a912465e9e
# ╟─4a50c199-2d06-420a-b797-3faa40ab14b8
# ╠═1e6e16d0-a4e4-43e0-a67b-88656c1f5247
# ╠═4664553e-d37a-41f7-a439-c833e8ba95a9
# ╠═8f3027d3-fb5f-4f24-ab6b-4dc4526a4f99
# ╠═f3203adb-6624-4cb8-bab7-825e9e2d956a
# ╠═79ded2f2-eb09-4114-abda-d79b6111fbb5
