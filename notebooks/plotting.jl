### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ bf726e2b-89fa-4cc0-8b12-bdac4e3d7580
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using PlutoUI
    using SimplexGridFactory
	using GridVisualize
	using ExtendableGrids
	using Triangulate
	using TetGen
	using PlutoVista
	default_plotter!(PlutoVista)
	TableOfContents()
end

# ╔═╡ df590d31-e543-4d3a-a3b1-4e7a411927a2
function condensator3D_tensorgridNEW(; scale = [50,50,50], d = 10, nrefs = 2)

    @info "Generating 3D condensator grid for a cuboid with dimensions
	($(scale[1]),$(scale[2]),$(2*scale[1]+d)) and middle layer width $d."

    X = scale[1]
    Y = scale[2]
    Z = 2*scale[3] + d
    npts = 2*(nrefs+1)

    XX = LinRange(0,scale[1],npts)
    YY = LinRange(0,scale[2],npts)

    xgrid = simplexgrid(XX,YY)

    z_levels_nonuniform = Array{Float64,1}(LinRange(0,scale[3]-d,npts))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3]-d,scale[3],1+Int(npts/2))[3:end]))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3],scale[3]+d,1+Int(npts/2))[2:end-1]))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3]+d,scale[3]+2*d,1+Int(npts/2))[1:end-2]))
    append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3]+2*d,Z,npts)))

    xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=4, top_offset=5)
    xgrid = uniform_refine(xgrid,nrefs-1)
    # the offsets lead to the following boundary regions:
    # 1 - 4  = side core
    # 5      = bottom core
    # 6      = tope core
    # 7 - 10 = side stressor

    cellmask!(xgrid,[0,0,scale[3]],[X,Y,scale[3]+d],2)
	bfacemask!(xgrid,[0,0,scale[3]],[X,0,scale[3]+d],7)
	bfacemask!(xgrid,[0,Y,scale[3]],[X,Y,scale[3]+d],8)
    bfacemask!(xgrid,[0,0,scale[3]],[0,Y,scale[3]+d],9)
	bfacemask!(xgrid,[X,0,scale[3]],[X,Y,scale[3]+d],10)

    return xgrid
end

# ╔═╡ 5644773c-c197-413c-9e31-b49bd8851996
begin
	xgrid_condensator3D_tensorgrid= condensator3D_tensorgridNEW(scale=[100,100,100], d=10, nrefs=2)
	gridplot(xgrid_condensator3D_tensorgrid; Plotter=PlutoVista)
end

# ╔═╡ d3b27958-25f0-43ec-bb57-7f07b30d7479
function condensator2D_tensorgrid(; scale = [50,50], d = 10, nrefs = 2)

    @info "Generating 2D condensator grid for a rectangle with dimensions
	($(scale[1]),$(2*scale[2]+d)) and middle layer width $d."

    X = scale[1]
    Y = 2*scale[2] + d
    npts = 2*(nrefs+1)

    XX = LinRange(0,scale[1],npts)
    YY = Array{Float64,1}(LinRange(0,scale[2]-d,npts))
    append!(YY, Array{Float64,1}(LinRange(scale[2]-d,scale[2],1+Int(npts/2))[3:end]))
    append!(YY, Array{Float64,1}(LinRange(scale[2],scale[2]+d,1+Int(npts/2))[2:end-1]))
    append!(YY, Array{Float64,1}(LinRange(scale[2]+d,scale[2]+2*d,1+Int(npts/2))[1:end-2]))
    append!(YY, Array{Float64,1}(LinRange(scale[2]+2*d,Y,npts)))

	xgrid = simplexgrid(XX,YY)
    xgrid = uniform_refine(xgrid,nrefs-1)
    # the offsets lead to the following boundary regions:
    # 1 = bottom bulk
    # 2 = right bulk
    # 3 = top bulk
    # 4 = left bulk
	# 5 = right side stressor
	# 6 = left side stressor

    cellmask!(xgrid,[0,scale[2]],[X,scale[2]+d],2)
    bfacemask!(xgrid,[X,scale[2]],[X,scale[2]+d],5)
	bfacemask!(xgrid,[0,scale[2]],[0,scale[2]+d],6)
	

    return xgrid
end

# ╔═╡ 9a54d847-5c35-417e-97ef-d7b62e60c2de
begin
	xgrid_condensator2D_tensorgrid= condensator2D_tensorgrid(scale=[100,100], d=10, nrefs=1)
	gridplot(xgrid_condensator2D_tensorgrid; Plotter=PlutoVista)
end

# ╔═╡ 4cf263cc-ab63-4358-b5c8-fc042d7d2968
function condensator2D_periodic(; A = 50, B = 100, d = 5, reflevel = 1, maxvol1 = B*d/4, maxvol2 = B*d/4)

    builder=SimplexGridBuilder(Generator=Triangulate)

    @info "Generating 2d condensator grid for A = $A, B = $B, d = $d"

    p1=point!(builder,0,0)                                                
    p2=point!(builder,B,0)                                         
    p3=point!(builder,B,A)                                 
    p4=point!(builder,B,A+d)    
    p5=point!(builder,B,2*A+d)                                                
    p6=point!(builder,0,2*A+d)                                         
    p7=point!(builder,0,A+d)                                 
    p8=point!(builder,0,A)     

    facetregion!(builder,1) 
    facet!(builder,p1,p2) # bottom of bottom plate
	
    facetregion!(builder,2)
    facet!(builder,p2,p3) # right boundary of bottom plate
    facet!(builder,p3,p4) # right gap boundary
    facet!(builder,p4,p5) # right boundary of top plate

	facetregion!(builder,3)
    facet!(builder,p5,p6) # top boundary of top plate
	
	facetregion!(builder,4)
    facet!(builder,p6,p7) # left boundary of top plate 
    facet!(builder,p7,p8) # left gap boundary
    facet!(builder,p8,p1) # left boundary of bottom plate

    # interior facets to split materials
	facetregion!(builder,6)
    facet!(builder,p3,p8)
    facet!(builder,p4,p7)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*B,0.5*A))
    regionpoint!(builder,(0.5*B,1.5*A+d))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*B,A+0.5*d))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end

# ╔═╡ abfb865e-e346-412e-ad3c-f44d25797bd2
begin
	xgrid_condensator2D = condensator2D_periodic(A=100,B=100,d = 20)
	gridplot(xgrid_condensator2D; Plotter=PlutoVista)
end

# ╔═╡ ff6cd0f3-3753-4c76-a249-007a2ad89b6e
function bimetal_tensorgrid_uniform(; scale = [1,1,1], nrefs = 1, material_border = 0.5)

    @info "Generating bimetal 3D grid for scale = $scale and middle interface at $material_border of height $scale[2]"

    W = scale[1]
    H = scale[2]
    Z = scale[3]
    h1 = round(scale[2]*material_border)
    h2 = round(scale[2]*(1 - material_border))
    hz_factor = 2.0^-nrefs

	hx = W/4*2.0^-nrefs
    hy = min(h1,h2)/2*2.0^-nrefs
	hz = 100 * hz_factor
    
	XX = 0:hx:W
	YY = Array{Float64,1}(0:(h1-hy)/2:h1-hy)
	append!(YY, Array{Float64,1}(LinRange(h1-hy,h1+hy,3)[2:end-1]))
	append!(YY, Array{Float64,1}(h1+hy:(h2-hy)/2:H))
	ZZ = Array{Float64,1}(0:hz:Z)

    xgrid = simplexgrid(XX,YY)
    xgrid = simplexgrid(xgrid, ZZ)
	xgrid = uniform_refine(xgrid,nrefs)
    # the offsets lead to the following boundary regions:
    # 1 - 6 = sides core & bottom
    # 7 = bottom core
    # 8 = bottom stressor
    # 9 = top core
	# 10 = top stressor

	# assigning region numbers: core region = 1, stressor region = 2
	cellmask!(xgrid,[0,0,0],[W,h1,Z],1)
	cellmask!(xgrid,[0,h1,0],[W,H,Z],2)

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

    return xgrid

end

# ╔═╡ c4f547f6-7064-4a0c-903b-bbd6b55bc140
begin
	xgrid_new = bimetal_tensorgrid_uniform(scale = [20,30,50], material_border = 0.7, nrefs = 1)
	gridplot(xgrid_new; Plotter=PlutoVista)
end

# ╔═╡ ff2ebfe6-ab76-470f-b96e-e74ebbbda2e1
function bimetal_tensorgrid(; scale = [1,1,1], nrefs = 1, material_border = 0.5)

    @info "Generating bimetal grid for scale = $scale"

    vol_factor = scale[1]*scale[2]*4.0^-nrefs/2
    hz_factor = 2.0^-nrefs

    builder=SimplexGridBuilder(Generator=Triangulate)
    p1=point!(builder,0,0)                                                
    p2=point!(builder,scale[1],0)                                         
    p3=point!(builder,scale[1],scale[2])                                  
    p4=point!(builder,0,scale[2])                                         
    p5=point!(builder,0,material_border*scale[2])                         
    p6=point!(builder,scale[1],material_border*scale[2])                 

    #facetregion!(builder,7) # interior boundary
    facet!(builder,p5 ,p6)

    facetregion!(builder,1) # outer boundary material 1
    facet!(builder,p5, p1)
    facet!(builder,p1, p2)
    facet!(builder,p2, p6)

    facetregion!(builder,2) # outer boundary material 2
    facet!(builder,p6, p3)
    facet!(builder,p3, p4)
    facet!(builder,p4, p5)

    cellregion!(builder,1)
    maxvolume!(builder,vol_factor)
    regionpoint!(builder,(0.5*scale[1],0.5*material_border*scale[2]))

    cellregion!(builder,2)
    maxvolume!(builder,vol_factor)
    regionpoint!(builder,(0.5*scale[1],0.5*(material_border+1)*scale[2]))

    xgrid = simplexgrid(builder)

    hz = 100 * hz_factor
    xgrid = simplexgrid(xgrid,0:hz:scale[3]; bot_offset = 2, top_offset = 4)

    return xgrid

end

# ╔═╡ 04b59757-fae4-4607-b36c-113ea85f0c57
begin
	xgrid_new2 = bimetal_tensorgrid(scale = [100,100,2000], material_border = 0.75, nrefs = 2)
	gridplot(xgrid_new2; Plotter=PlutoVista)
end

# ╔═╡ be4abc28-126f-4578-aa81-674d7d096403
function condensator3D(; scale=[50,50,50], d = 5, nrefs = 1, maxvol1 = prod(scale)/4, maxvol2 = scale[1]*scale[2]*d/4)

    builder=SimplexGridBuilder(Generator=TetGen)

    @info "Generating 3D condensator grid for scale = $scale and middle layer width = $d."

	X = scale[1]
	Y = scale[2]
	Z = scale[3]
	
    # side at y = 0
    p1=point!(builder,0,0,0)
    p2=point!(builder,X,0,0)
    p3=point!(builder,X,0,Z)
    p4=point!(builder,X,0,Z+d)
    p5=point!(builder,X,0,2*Z+d)
    p6=point!(builder,0,0,2*Z+d)
    p7=point!(builder,0,0,Z+d)
    p8=point!(builder,0,0,Z)

    # side at y = Y
    p9=point!(builder,0,Y,0)
    p10=point!(builder,X,Y,0)
    p11=point!(builder,X,Y,Z)
    p12=point!(builder,X,Y,Z+d)
    p13=point!(builder,X,Y,2*Z+d)
    p14=point!(builder,0,Y,2*Z+d)
    p15=point!(builder,0,Y,Z+d)
    p16=point!(builder,0,Y,Z)

    facetregion!(builder,1) # bottom of bottom cube
    facet!(builder,p1,p9,p10,p2)
    facetregion!(builder,2) # front, back, left and right boundary of bottom cube
    facet!(builder,p1,p2,p3,p8)
    facet!(builder,p10,p9,p16,p11)
    facet!(builder,p2,p10,p11,p3)
    facet!(builder,p9,p1,p8,p16)

    facetregion!(builder,3) # gap boundary front, back
    facet!(builder,p8,p3,p4,p7)
    facet!(builder,p11,p16,p15,p12)
    facetregion!(builder,4) # gap boundary left and right
    facet!(builder,p3,p11,p12,p4)
    facet!(builder,p16,p8,p7,p15)

    facetregion!(builder,5) # boundary of top cube (front, back, left, right)
    facet!(builder,p7,p4,p5,p6)
    facet!(builder,p12,p15,p14,p13)
    facet!(builder,p4,p12,p13,p5)
    facet!(builder,p15,p7,p6,p14)
    facetregion!(builder,6) # boundary of top cube (top)
    facet!(builder,p6,p5,p13,p14)

    #facetregion!(builder,7) # interior facets to split materials
    facet!(builder,p8,p3,p11,p16)
    facet!(builder,p7,p4,p12,p15)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*X,0.5*Y,0.5*Z))
    regionpoint!(builder,(0.5*X,0.5*Y,1.5*Z+d))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*X,0.5*Y,Z+0.5*d))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,nrefs)

    return xgrid
end

# ╔═╡ 9e18c4db-0db3-4a1d-8eb5-f9668f8e7f7e
function condensator3D_tensorgrid_new(; scale = [50,50,50], d = 10, nrefs = 2)
	
	@info "Generating 3D condensator tensor grid for scale = $scale and middle layer width = $d."

	X = scale[1]
    Y = scale[2]
    Z = 2*scale[3] + d
    hz_factor = 2.0^-nrefs

	hx = 100*2.0^-nrefs
	hy = hx
	XX = 0:hx:X
	YY = 0:hy:Y

	xgrid = simplexgrid(XX,YY)

	hz = 50 * hz_factor
	z_levels_nonuniform = Array{Float64,1}(0:hz:scale[3])
	if nrefs == 1
		npts = 2*nrefs
	else
		npts = nrefs
	end
	append!(z_levels_nonuniform, Array{Float64,1}(LinRange(scale[3],scale[3]+d,npts)[2:end-1]))
	append!(z_levels_nonuniform, Array{Float64,1}(scale[3]+d:hz:Z))

    xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=4, top_offset=5)
	xgrid = uniform_refine(xgrid,nrefs-1)
    # the offsets lead to the following boundary regions:
	# 1 - 4 = side core
	# 5 = bottom core
	# 6 = tope core
	# 7 = side stressor

	cellmask!(xgrid,[0.0,0.0,scale[3]],[X,Y,scale[3]+d],2)
    bfacemask!(xgrid,[0.0,0.0,scale[3]],[X,Y,scale[3]+d],7)
	
    return xgrid
end

# ╔═╡ 272e025b-c6ee-43a2-99a3-5a88211427c2
begin
	xgrid_1 = condensator3D_tensorgrid_new(; scale = [500,500,500], d = 10, nrefs = 1)
	gridplot(xgrid_1; Plotter=PlutoVista)
end

# ╔═╡ ed12f497-aee9-4003-b110-bc032d5b4725
begin
	xgrid1 = condensator3D(; scale = [100,100,100], d = 10, nrefs = 2)
	gridplot(xgrid1; Plotter=PlutoVista)
end

# ╔═╡ 0f48688f-e4b5-4628-b90a-9a3c3c4266d0
function condensator3D_tensorgrid(; scale = [50,50,50], d = 10, nrefs = 1)
	
	@info "Generating 3D condensator tensor grid for scale = $scale and middle layer width = $d."

	builder=SimplexGridBuilder(Generator=Triangulate)

	X = scale[1]
    Y = scale[2]
    Z = scale[3]
    vol_factor_core = 4.0^-nrefs
    vol_factor_stressor = 4.0^-nrefs
    hz_factor = 2.0^-nrefs
	maxvol1 = X*Y*vol_factor_core
	maxvol2 = X*d*vol_factor_core

    # side at y = 0
	#p0 = point!(builder,0.5*X,0.5*Y)
    p1 = point!(builder,0,0)
    p2 = point!(builder,X,0)
    #p3 = point!(builder,X,Y-2*d)
	p4 = point!(builder,X,Y)
    p5 = point!(builder,X,Y+d)
	#p6 = point!(builder,X,Y+3*d)
    p7 = point!(builder,X,2*Y+d)
    p8 = point!(builder,0,2*Y+d)
	#p9 = point!(builder,0,Y+3*d)
    p10 = point!(builder,0,Y+d)
    p11 = point!(builder,0,Y)
	#p12 = point!(builder,0,Y-2*d)
	#p13 = point!(builder,0.5*X,1.5*Y+d)
	p14 = point!(builder,0.5*X,Y)
	p15 = point!(builder,0.5*X,Y+d)

	facetregion!(builder,1) # bottom rectangle
	#facet!(builder,p0,p1)
	#facet!(builder,p0,p2)
	#facet!(builder,p0,p3)
	#facet!(builder,p0,p12)
	facet!(builder,p1,p2)
	facet!(builder,p2,p4)
	facet!(builder,p4,p11)
	facet!(builder,p11,p1)

	facetregion!(builder,2) # side of middle layer
    facet!(builder,p4,p5)
	facet!(builder,p10,p11)
	
	facetregion!(builder,1) # top rectangle
	#facet!(builder,p9,p7)
	#facet!(builder,p9,p4)
	#facet!(builder,p9,p5)
	#facet!(builder,p9,p6)
    facet!(builder,p5,p7)
	facet!(builder,p7,p8)
	facet!(builder,p8,p10)
	facet!(builder,p10,p5)

	cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
	regionpoint!(builder,(0.5*X,0.5*Y))
	regionpoint!(builder,(0.5*X,1.5*Y+d))

    cellregion!(builder,2) # material 2
	maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*X,Y+0.5*d))

	xgrid = simplexgrid(builder)

    hz = 100/nrefs * hz_factor
    z_levels = 0:hz:Z
    z_levels_nonuniform = Vector{Any}(z_levels)
    z_levels_nonuniform = vcat(z_levels_nonuniform...)

    xgrid = simplexgrid(xgrid, z_levels_nonuniform; bot_offset=2, top_offset=4)
	#xgrid = uniform_refine(xgrid,nrefs)
    # the offsets lead to the following boundary regions:
    # 1 = side core (not seen from outside)
    # 2 = side stressor
    # 3 = bottom core
    # 4 = bottom stressor
    # 5 = top core
    # 6 = top stressor
	
    return xgrid
end

# ╔═╡ 49dde1de-dcfe-47df-8a35-e53ea25c366b
begin
	xgridnew = condensator3D_tensorgrid(; scale=[100,100,100], d = 10, nrefs = 2)
	gridplot(xgridnew; Plotter=PlutoVista)
end

# ╔═╡ 3e5f4ebc-38e7-4067-a0c7-478829d3f1b2
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

# ╔═╡ 634c46b2-1a15-46c3-ad8d-57d9b0ea5695
geometry = [30, 20, 10, 2000]

# ╔═╡ f0e17b2a-38d4-475f-ada7-3277542d64b6
begin
	xgridnew3 = nanowire_tensorgrid(; scale = geometry, cut_levels = geometry[4]/2, nrefs = 2, α = nothing, version = 1)
	
	gridplot(xgridnew3; Plotter=PlutoVista)
end

# ╔═╡ c024ebaf-f408-4cb8-81e1-8dca6cd8262c
begin
	xgrid = nanowire_tensorgrid(; scale = geometry, cut_levels = geometry[4]/2, nrefs = 2, α = geometry[3]/4, version = 1)
	
	gridplot(xgrid; Plotter=PlutoVista)
end

# ╔═╡ 22e14399-c68f-4fc4-b55f-35cf2f052725
begin
	xgrid2 = nanowire_tensorgrid(; scale = geometry, cut_levels = geometry[4]/2, nrefs = 1, α = nothing)
	
	gridplot(xgrid2; Plotter=PlutoVista)
end

# ╔═╡ 9a4b51d6-1f3c-4777-b094-546d71ee2c46
begin
	xgrid3 = nanowire_tensorgrid(; scale = [25/sqrt(3),20/sqrt(3),20,500], cut_levels = nothing, nrefs = 0, α = nothing,z_levels_dist = 100, version = 1)
	gridplot(xgrid3; Plotter=PlutoVista)
end

# ╔═╡ 4367a612-8853-4f7b-ad0f-240741525215
begin
	xgrid4 = nanowire_tensorgrid(; scale = [20/sqrt(3),10/sqrt(3),10,500], cut_levels = nothing, nrefs = 2, α = nothing,z_levels_dist = 100, version = 2)
	gridplot(xgrid4; Plotter=PlutoVista)
end

# ╔═╡ 70e49981-53de-4e91-b38f-5c67368a2370
begin
	xgrid5 = nanowire_tensorgrid(; scale = [20/sqrt(3),30/sqrt(3),30,1000], cut_levels = nothing, nrefs = 1, α = nothing,z_levels_dist = 25, version = 2)
	gridplot(xgrid5; Plotter=PlutoVista)
end

# ╔═╡ Cell order:
# ╠═bf726e2b-89fa-4cc0-8b12-bdac4e3d7580
# ╠═df590d31-e543-4d3a-a3b1-4e7a411927a2
# ╠═5644773c-c197-413c-9e31-b49bd8851996
# ╠═d3b27958-25f0-43ec-bb57-7f07b30d7479
# ╠═9a54d847-5c35-417e-97ef-d7b62e60c2de
# ╠═4cf263cc-ab63-4358-b5c8-fc042d7d2968
# ╠═abfb865e-e346-412e-ad3c-f44d25797bd2
# ╠═ff6cd0f3-3753-4c76-a249-007a2ad89b6e
# ╠═c4f547f6-7064-4a0c-903b-bbd6b55bc140
# ╠═ff2ebfe6-ab76-470f-b96e-e74ebbbda2e1
# ╠═04b59757-fae4-4607-b36c-113ea85f0c57
# ╠═be4abc28-126f-4578-aa81-674d7d096403
# ╠═9e18c4db-0db3-4a1d-8eb5-f9668f8e7f7e
# ╠═272e025b-c6ee-43a2-99a3-5a88211427c2
# ╠═ed12f497-aee9-4003-b110-bc032d5b4725
# ╠═0f48688f-e4b5-4628-b90a-9a3c3c4266d0
# ╠═49dde1de-dcfe-47df-8a35-e53ea25c366b
# ╠═3e5f4ebc-38e7-4067-a0c7-478829d3f1b2
# ╠═634c46b2-1a15-46c3-ad8d-57d9b0ea5695
# ╠═f0e17b2a-38d4-475f-ada7-3277542d64b6
# ╠═c024ebaf-f408-4cb8-81e1-8dca6cd8262c
# ╠═22e14399-c68f-4fc4-b55f-35cf2f052725
# ╠═9a4b51d6-1f3c-4777-b094-546d71ee2c46
# ╠═4367a612-8853-4f7b-ad0f-240741525215
# ╠═70e49981-53de-4e91-b38f-5c67368a2370
