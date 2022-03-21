function bimetal_strip3D(; material_border = 0.5, scale = [1,1,1], anisotropy = [1,1,2], reflevel = 1, maxvol = prod(scale./anisotropy)/8)

    @info "Generating bimetal grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=TetGen)
    p1=point!(builder,0,0,0)                                                
    p2=point!(builder,scale[1],0,0)                                         
    p3=point!(builder,scale[1],scale[2],0)                                  
    p4=point!(builder,0,scale[2],0)                                         
    p5=point!(builder,0,0,scale[3])                                         
    p6=point!(builder,scale[1],0,scale[3])                                  
    p7=point!(builder,scale[1],scale[2],scale[3])                           
    p8=point!(builder,0,scale[2],scale[3])                                  

    p9=point!(builder,0,material_border*scale[2],0)                         
    p10=point!(builder,scale[1],material_border*scale[2],0)                 
    p11=point!(builder,scale[1],material_border*scale[2],scale[3])          
    p12=point!(builder,0,material_border*scale[2],scale[3])                 

    facetregion!(builder,1) # bottom (material A) = core front
    facet!(builder,p1 ,p2 ,p10 ,p9)
    facet!(builder,p9 ,p10 ,p3 ,p4)

    facetregion!(builder,2) # shell boundary
    facet!(builder,p6 ,p2 ,p1 ,p5)
    facet!(builder,p11 ,p6 ,p5 ,p12)
    facet!(builder,p2 ,p6 ,p11 ,p10)
    facet!(builder,p1 ,p9 ,p12 ,p5)

    facetregion!(builder,3) # stressor boundary
    facet!(builder,p7 ,p11 ,p12 ,p8)
    facet!(builder,p10 ,p11 ,p7 ,p3)
    facet!(builder,p3 ,p7 ,p8 ,p4)
    facet!(builder,p9 ,p4 ,p8 ,p12)

    facetregion!(builder,4) # interior facet to split regions
    facet!(builder,p9 ,p10 ,p11 ,p12)

    cellregion!(builder,1)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[1],0.5*material_border*scale[2],0.5*scale[3]))

    cellregion!(builder,2)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[1],0.5*(material_border+1)*scale[2],0.5*scale[3]))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 3, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)
    return xgrid

end


function bimetal_strip3D_middle_layer(; material_border = 0.5, scale = [1,1,1], anisotropy = [1,1,2], reflevel = 1, maxvol1 = prod(scale./anisotropy)/8, maxvol2 = prod(scale./anisotropy))

    @info "Generating bimetal grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=TetGen)
    y = 5
    p1=point!(builder,0,0,0)
    p2=point!(builder,scale[1],0,0)
    p3=point!(builder,scale[1],material_border*scale[2]-y,0)
    p4=point!(builder,scale[1],material_border*scale[2],0)
    p5=point!(builder,scale[1],material_border*scale[2]+y,0)
    p6=point!(builder,scale[1],scale[2],0)
    p7=point!(builder,0,scale[2],0)
    p8=point!(builder,0,material_border*scale[2]+y,0)
    p9=point!(builder,0,material_border*scale[2],0)
    p10=point!(builder,0,material_border*scale[2]-y,0)

    p11=point!(builder,0,0,scale[3])
    p12=point!(builder,scale[1],0,scale[3])
    p13=point!(builder,scale[1],material_border*scale[2]-y,scale[3])
    p14=point!(builder,scale[1],material_border*scale[2],scale[3])
    p15=point!(builder,scale[1],material_border*scale[2]+y,scale[3])
    p16=point!(builder,scale[1],scale[2],scale[3])
    p17=point!(builder,0,scale[2],scale[3])
    p18=point!(builder,0,material_border*scale[2]+y,scale[3])
    p19=point!(builder,0,material_border*scale[2],scale[3])
    p20=point!(builder,0,material_border*scale[2]-y,scale[3])

    facetregion!(builder,1) # core boundary (bottom, top, sides)
    facet!(builder,p1,p2,p3,p10)
    facetregion!(builder,12)
    facet!(builder,p11,p12,p13,p20)
    facetregion!(builder,13)
    facet!(builder,p1,p2,p12,p11)
    facet!(builder,p2,p3,p13,p12)
    facet!(builder,p1,p10,p20,p11)

    facetregion!(builder,2) # middle layer 1 (bottom, top, sides)
    facet!(builder,p3,p4,p9,p10)
    facet!(builder,p13,p14,p19,p20)
    facet!(builder,p3,p10,p20,p13)
    facet!(builder,p3,p4,p14,p13)
    facet!(builder,p9,p10,p20,p19)

    facetregion!(builder,3) # middle layer 2 (bottom, top, sides)
    facet!(builder,p4,p5,p8,p9)
    facet!(builder,p14,p15,p18,p19)
    facet!(builder,p5,p8,p18,p15)
    facet!(builder,p4,p5,p15,p14)
    facet!(builder,p8,p9,p19,p18)

    facetregion!(builder,4) # stressor boundary (bottom, top, sides)
    facet!(builder,p5,p6,p7,p8)
    facet!(builder,p15,p16,p17,p18)
    facet!(builder,p6,p7,p17,p16)
    facet!(builder,p5,p6,p16,p15)
    facet!(builder,p7,p8,p18,p17)

    facetregion!(builder,5) # interior facet to split regions
    facet!(builder,p4,p9,p19,p14)

    cellregion!(builder,1)
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*scale[1],0.5*(material_border*scale[2]-y),0.5*scale[3]))
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*scale[1],material_border*scale[2]-y/2,0.5*scale[3]))

    cellregion!(builder,2)
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*scale[1],0.5*((material_border+1)*scale[2]+y),0.5*scale[3]))
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*scale[1],material_border*scale[2]+y/2,0.5*scale[3]))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 3, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)
    return xgrid

end


function bimetal_strip2D(; material_border = 0.5, scale = [1,1], anisotropy = [1,1], reflevel = 1, maxvol = prod(scale./anisotropy)/4)

    @info "Generating 2d bimetal grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=Triangulate)
    p1=point!(builder,0,0)                                                
    p2=point!(builder,scale[2],0)                                         
    p3=point!(builder,scale[2],scale[1])                                 
    p4=point!(builder,0,scale[1])                                      

    p5=point!(builder,0,material_border*scale[1])                         
    p6=point!(builder,scale[2],material_border*scale[1])


    facetregion!(builder,1) # left (material A)
    facet!(builder,p5 ,p1)
    facetregion!(builder,11) # left (material B)
    facet!(builder,p4,p5)

    facetregion!(builder,2) # bottom
    facet!(builder,p1 ,p2)

    facetregion!(builder,3) # top
    facet!(builder,p3, p4)

    facetregion!(builder,4) # right
    facet!(builder,p2, p6)
    facet!(builder,p6, p3)

    facetregion!(builder,99) # interior facet to split regions
    facet!(builder,p5, p6)

    cellregion!(builder,1)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[2],0.5*material_border*scale[1]))

    cellregion!(builder,2)
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0.5*scale[2],0.5*(material_border+1)*scale[1]))

    xgrid = simplexgrid(builder)

    Coords = xgrid[Coordinates]
    for k = 1 : 2, n = 1 : size(Coords,2)
        Coords[k,n] *= anisotropy[k]
    end
    scale .*= anisotropy

    xgrid = uniform_refine(xgrid,reflevel)
    return xgrid
    
end


function condensator3D(; A = 50, B = 100, Z = 20, d = 5, reflevel = 1, maxvol1 = Z*B*d/4, maxvol2 = Z*B*d/4)

    builder=SimplexGridBuilder(Generator=TetGen)

    @info "Generating 3d condensator grid for A = $A, B = $B, Z = $Z, d = $d"

    # bottom side at Z = 0
    p1=point!(builder,0,0,0)
    p2=point!(builder,B,0,0)
    p3=point!(builder,B,0,A)
    p4=point!(builder,B,0,A+d)
    p5=point!(builder,B,0,2*A+d)
    p6=point!(builder,0,0,2*A+d)
    p7=point!(builder,0,0,A+d)
    p8=point!(builder,0,0,A)

    # top side at Z
    p9=point!(builder,0,Z,0)
    p10=point!(builder,B,Z,0)
    p11=point!(builder,B,Z,A)
    p12=point!(builder,B,Z,A+d)
    p13=point!(builder,B,Z,2*A+d)
    p14=point!(builder,0,Z,2*A+d)
    p15=point!(builder,0,Z,A+d)
    p16=point!(builder,0,Z,A)

    facetregion!(builder,1) # bottom of bottom cube
    facet!(builder,p1 ,p9, p10, p2)
    facetregion!(builder,2) # front, back, left and right boundary of bottom cube
    facet!(builder,p1,p2,p3,p8)
    facet!(builder,p10,p9,p16,p11)
    facet!(builder,p2,p10,p11,p3)
    facet!(builder,p9,p1,p8,p16)

    facetregion!(builder,2) # gap boundary front, back, left and right
    facet!(builder,p8 ,p3, p4, p7)
    facet!(builder,p11 ,p16, p15, p12)
    facet!(builder,p3 ,p11, p12, p4)
    facet!(builder,p16 ,p8, p7, p15)

    facetregion!(builder,3) # boundary of top cube (front, back, left, top, right)
    facet!(builder,p7, p4, p5, p6)
    facet!(builder,p12, p15, p14, p13)
    facet!(builder,p4, p12, p13, p5)
    facet!(builder,p6, p5, p13, p14)
    facet!(builder,p15,p7,p6,p14)

    facetregion!(builder,99) # interior facets to split materials
    facet!(builder,p8, p3, p11, p16)
    facet!(builder,p7, p4, p12, p15)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*B,0.5*Z,0.5*A))
    regionpoint!(builder,(0.5*B,0.5*Z,1.5*A+d))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,(0.5*B,0.5*Z,A+0.5*d))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function condensator2D(; A = 50, B = 100, d = 5, reflevel = 1, maxvol1 = B*d/4, maxvol2 = B*d/4)

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

    facetregion!(builder,1) # bottom of bottom plate
    facet!(builder,p1 ,p2)
    facetregion!(builder,11) # left and right boundary of bottom plate
    facet!(builder,p2,p3)
    facet!(builder,p8,p1)

    facetregion!(builder,2) # gap boundary left and right
    facet!(builder,p3 ,p4)
    facet!(builder,p7 ,p8)

    facetregion!(builder,3) # boundary of top plate (left, top, right)
    facet!(builder,p4, p5)
    facet!(builder,p5, p6)
    facet!(builder,p6, p7)

    facetregion!(builder,4) # interior facets to split materials
    facet!(builder,p3, p8)
    facet!(builder,p4, p7)

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

function nonpolarquantumwell3D(; X = 80, Y = 50, Z = 20, V = 8, d = 2, w = 3, reflevel = 1, maxvol1 = 1/4*(Y*Z*(X-V)/2), maxvol2 = 1/4*(Y*Z*V))

    builder=SimplexGridBuilder(Generator=TetGen)

    @info "Generating 3d quantum well grid for X = $X, Y = $Y, Z = $Z, V = $V, d = $d, w = $w"

    # bottom side at Y = 0
    p1=point!(builder,0,0,0)
    p2=point!(builder,0,0,Z)
    p3=point!(builder,(X-V)/2,0,Z)
    p4=point!(builder,(X-V)/2+V,0,Z)
    p5=point!(builder,X,0,Z)
    p6=point!(builder,X,0,0)
    p7=point!(builder,(X-V)/2+V,0,0)
    p8=point!(builder,(X-V)/2,0,0)
    p9=point!(builder,(X-V)/2+V,0,(Z+w)/2)
    p10=point!(builder,(X-V)/2+V-d,0,(Z+w)/2)
    p11=point!(builder,(X-V)/2+V-d,0,(Z-w)/2)
    p12=point!(builder,(X-V)/2+V,0,(Z-w)/2)

    # top side at Y
    p13=point!(builder,0,Y,0)
    p14=point!(builder,0,Y,Z)
    p15=point!(builder,(X-V)/2,Y,Z)
    p16=point!(builder,(X-V)/2+V,Y,Z)
    p17=point!(builder,X,Y,Z)
    p18=point!(builder,X,Y,0)
    p19=point!(builder,(X-V)/2+V,Y,0)
    p20=point!(builder,(X-V)/2,Y,0)
    p21=point!(builder,(X-V)/2+V,Y,(Z+w)/2)
    p22=point!(builder,(X-V)/2+V-d,Y,(Z+w)/2)
    p23=point!(builder,(X-V)/2+V-d,Y,(Z-w)/2)
    p24=point!(builder,(X-V)/2+V,Y,(Z-w)/2)

    facetregion!(builder,1) # sides of region 1 cube (bottom, top, front, left, right)
    facet!(builder,p1,p2,p3,p8)
    facet!(builder,p13,p14,p15,p20)
    facet!(builder,p1,p2,p14,p13)
    facet!(builder,p1,p8,p20,p13)
    facet!(builder,p2,p3,p15,p14)

    facetregion!(builder,2) # gap region
    # bottom
    facet!(builder,p3,p4,p9,p10)
    facet!(builder,p3,p10,p11,p8)
    facet!(builder,p8,p11,p12,p7)
    # top
    facet!(builder,p15,p16,p21,p22)
    facet!(builder,p15,p22,p23,p20)
    facet!(builder,p20,p23,p24,p19)
    # front
    facet!(builder,p3,p8,p20,p15)
    # back
    facet!(builder,p4,p9,p21,p16)
    facet!(builder,p9,p10,p22,p21)
    facet!(builder,p10,p11,p23,p22)
    facet!(builder,p11,p12,p24,p23)
    facet!(builder,p12,p7,p19,p24)
    # left
    facet!(builder,p7,p8,p20,p19)
    # right
    facet!(builder,p3,p4,p16,p15)

    facetregion!(builder,3) # sides of region 2 cube
    # bottom
    facet!(builder,p4,p5,p9)
    facet!(builder,p5,p9,p12,p6)
    facet!(builder,p6,p12,p7)
    #facet!(builder,p4,p5,p6,p9)
    #facet!(builder,p6,p9,p12,p7)
    facet!(builder,p9,p10,p11,p12)
    # top
    facet!(builder,p16,p17,p18,p21)
    facet!(builder,p18,p21,p24,p19)
    facet!(builder,p21,p22,p23,p24)
    # back
    facet!(builder,p5,p6,p18,p17)
    # left
    facet!(builder,p6,p7,p19,p18)
    # right
    facet!(builder,p4,p5,p17,p16)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*(X-V)/2,0.5*Y,0.5*Z))
    regionpoint!(builder,(X-0.5*(X-V)/2,0.5*Y,0.5*Z))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,((X-V)/2+0.5*d,0.5*Y,0.5*Z))

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function nonpolarquantumwell2D(; X = 80, Z = 20, V = 8, d = 2, w = 3, reflevel = 1, maxvol1 = 1/4*(Z*(X-V)/2), maxvol2 = 1/4*(Z*V))

    builder=SimplexGridBuilder(Generator=Triangulate)

    @info "Generating 2d quantum well grid for X = $X, Z = $Z, V = $V, d = $d, w = $w"

    p1=point!(builder,0,0)
    p2=point!(builder,0,Z)
    p3=point!(builder,(X-V)/2,Z)
    p4=point!(builder,(X-V)/2+V,Z)
    p5=point!(builder,X,Z)
    p6=point!(builder,X,0)
    p7=point!(builder,(X-V)/2+V,0)
    p8=point!(builder,(X-V)/2,0)
    p9=point!(builder,(X-V)/2+V,(Z+w)/2)
    p10=point!(builder,(X-V)/2+V-d,(Z+w)/2)
    p11=point!(builder,(X-V)/2+V-d,(Z-w)/2)
    p12=point!(builder,(X-V)/2+V,(Z-w)/2)

    facetregion!(builder,1) # boundary of bottom plate (bottom, right, left)
    facet!(builder,p1,p2)
    facet!(builder,p2,p3)
    facet!(builder,p8,p1)

    facetregion!(builder,2) # gap boundary (right, left)
    facet!(builder,p3,p4)
    facet!(builder,p7,p8)

    facetregion!(builder,3) # interior facets to split materials
    facet!(builder,p3,p8)
    facet!(builder,p4,p9)
    facet!(builder,p9,p10)
    facet!(builder,p10,p11)
    facet!(builder,p11,p12)
    facet!(builder,p12,p7)

    facetregion!(builder,4) # boundary of top plate (left, top, right)
    facet!(builder,p4,p5)
    facet!(builder,p5,p6)
    facet!(builder,p6,p7)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol1)
    regionpoint!(builder,(0.5*(X-V)/2),0.5*Z)
    regionpoint!(builder,(X-0.5*(X-V)/2),0.5*Z)

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol2)
    regionpoint!(builder,((X-V)/2+0.5*d),0.5*Z)

    xgrid = simplexgrid(builder)
    xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end


function nanowire_grid(; scale = [1,1,1,1], anisotropy = [1,1,1,1], reflevel = 1, maxvol = prod(scale./anisotropy)/8)

    @info "Generating nanowire grid for scale = $scale"
    scale ./= anisotropy
    builder=SimplexGridBuilder(Generator=TetGen)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ = scale[3]

    # bottom side at Z = 0
    p0=point!(builder,0,0,0)
    p1=point!(builder,d1,0,0)
    p2=point!(builder,d1/2,sqrt(3)/2*d1,0)
    p3=point!(builder,-d1/2,sqrt(3)/2*d1,0)
    p4=point!(builder,-d1,0,0)
    p5=point!(builder,-d1/2,-sqrt(3)/2*d1,0)
    p6=point!(builder,d1/2,-sqrt(3)/2*d1,0)
    p7=point!(builder,d2,0,0)
    p8=point!(builder,d2/2,sqrt(3)/2*d2,0)
    p9=point!(builder,-d2/2,sqrt(3)/2*d2,0)
    p10=point!(builder,-d2,0,0)
    p11=point!(builder,-d2/2,-sqrt(3)/2*d2,0)
    p12=point!(builder,d2/2,-sqrt(3)/2*d2,0)
    p13=point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ,0)
    p14=point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ,0)

    # top side at Z = scale[4]
    p15=point!(builder,d1,0,scale[4])
    p16=point!(builder,d1/2,sqrt(3)/2*d1,scale[4])
    p17=point!(builder,-d1/2,sqrt(3)/2*d1,scale[4])
    p18=point!(builder,-d1,0,scale[4])
    p19=point!(builder,-d1/2,-sqrt(3)/2*d1,scale[4])
    p20=point!(builder,d1/2,-sqrt(3)/2*d1,scale[4])
    p21=point!(builder,d2,0,scale[4])
    p22=point!(builder,d2/2,sqrt(3)/2*d2,scale[4])
    p23=point!(builder,-d2/2,sqrt(3)/2*d2,scale[4])
    p24=point!(builder,-d2,0,scale[4])
    p25=point!(builder,-d2/2,-sqrt(3)/2*d2,scale[4])
    p26=point!(builder,d2/2,-sqrt(3)/2*d2,scale[4])
    p27=point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ,scale[4])
    p28=point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ,scale[4])
    p29=point!(builder,0,0,scale[4])

    facetregion!(builder,1) # core region (bottom, sides, top)
    # bottom
    facet!(builder,p0,p1,p2)
    facet!(builder,p0,p2,p3)
    facet!(builder,p0,p3,p4)
    facet!(builder,p0,p4,p5)
    facet!(builder,p0,p5,p6)
    facet!(builder,p0,p6,p1)
    #facet!(builder,p1,p2,p3,p4)
    #facet!(builder,p4,p5,p6,p1)
    # sides
    facetregion!(builder,11) # core region (bottom, sides, top)
    facet!(builder,p1,p2,p16,p15)
    facet!(builder,p2,p3,p17,p16)
    facet!(builder,p3,p4,p18,p17)
    facet!(builder,p4,p5,p19,p18)
    facet!(builder,p5,p6,p20,p19)
    facet!(builder,p6,p1,p15,p20)
    # top
    facetregion!(builder,12) # core region (bottom, sides, top)
    facet!(builder,p29,p15,p16)
    facet!(builder,p29,p16,p17)
    facet!(builder,p29,p17,p18)
    facet!(builder,p29,p18,p19)
    facet!(builder,p29,p19,p20)
    facet!(builder,p29,p20,p15)
    #facet!(builder,p15,p16,p17,p18)
    #facet!(builder,p18,p19,p20,p15)

    facetregion!(builder,2) # shell region (bottom, sides, top)
    # bottom
    facet!(builder,p1,p7,p8,p2)
    facet!(builder,p2,p8,p9,p3)
    facet!(builder,p3,p9,p10,p4)
    facet!(builder,p4,p10,p11,p5)
    facet!(builder,p5,p11,p12,p6)
    facet!(builder,p6,p12,p7,p1)
    # sides
    facet!(builder,p7,p8,p22,p21)
    facet!(builder,p8,p9,p23,p22)
    facet!(builder,p9,p10,p24,p23)
    facet!(builder,p10,p11,p25,p24)
    facet!(builder,p11,p12,p26,p25)
    facet!(builder,p12,p7,p21,p26)
    # top
    facet!(builder,p15,p21,p22,p16)
    facet!(builder,p16,p22,p23,p17)
    facet!(builder,p17,p23,p24,p18)
    facet!(builder,p18,p24,p25,p19)
    facet!(builder,p19,p25,p26,p20)
    facet!(builder,p20,p26,p21,p15)

    facetregion!(builder,3) # stressor region (bottom, sides, top)
    # bottom
    facet!(builder,p7,p12,p13)
    facet!(builder,p11,p12,p13,p14)
    facet!(builder,p10,p11,p14)
    # sides
    facet!(builder,p7,p13,p27,p21)
    facet!(builder,p13,p14,p28,p27)
    facet!(builder,p10,p14,p28,p24)
    # top
    facet!(builder,p21,p26,p27)
    facet!(builder,p25,p26,p27,p28)
    facet!(builder,p24,p25,p28)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0,0,scale[4]/2))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(scale[1]+scale[2]/2,0,scale[4]/2))

    cellregion!(builder,3) # material 3
    maxvolume!(builder,maxvol)
    regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2,scale[4]/2))

    xgrid = simplexgrid(builder)

    # Coords = xgrid[Coordinates]
    # for k = 1 : 3, n = 1 : size(Coords,2)
    #     Coords[k,n] *= anisotropy[k]
    # end
    # scale .*= anisotropy

    # xgrid = uniform_refine(xgrid,reflevel)

    return xgrid
end

function nanowire_tensorgrid(; scale = [1,1,1,1], nrefs = 1)

    @info "Generating nanowire grid for scale = $scale"
    
    builder=SimplexGridBuilder(Generator=Triangulate)

    d1 = scale[1]
    d2 = scale[1] + scale[2]
    δ = scale[3]
    vol_factor_core = 4.0^-nrefs
    vol_factor_shell = 4.0^-nrefs
    vol_factor_stressor = 4.0^-(nrefs+1)
    hz_factor = 2.0^-nrefs

    A_core = 3*sqrt(3)/2 * scale[1]^2
    A_shell = 3*sqrt(3)/2 * (scale[2]^2 + 2*scale[1]*scale[2])
    A_stressor = sqrt(3)/2 * scale[3] * (7*(scale[1]+scale[2]) + 3*scale[3])

    # bottom side at Z = 0
    #p0=point!(builder,0,0)
    p1=point!(builder,d1,0)
    p2=point!(builder,d1/2,sqrt(3)/2*d1)
    p3=point!(builder,-d1/2,sqrt(3)/2*d1)
    p4=point!(builder,-d1,0)
    p5=point!(builder,-d1/2,-sqrt(3)/2*d1)
    p6=point!(builder,d1/2,-sqrt(3)/2*d1)
    p7=point!(builder,d2,0)
    p8=point!(builder,d2/2,sqrt(3)/2*d2)
    p9=point!(builder,-d2/2,sqrt(3)/2*d2)
    p10=point!(builder,-d2,0)
    p11=point!(builder,-d2/2,-sqrt(3)/2*d2)
    p12=point!(builder,d2/2,-sqrt(3)/2*d2)
    p13=point!(builder,d2/2+δ/sqrt(3),-sqrt(3)/2*d2-δ)
    p14=point!(builder,-d2/2-δ/sqrt(3),-sqrt(3)/2*d2-δ)

    facetregion!(builder,1) # core region (bottom, sides, top)
    # bottom
    facet!(builder,p1,p2)
    facet!(builder,p2,p3)
    facet!(builder,p3,p4)
    facet!(builder,p4,p5)
    facet!(builder,p5,p6)
    facet!(builder,p6,p1)

    facetregion!(builder,2) # shell region (bottom, sides, top)
    # bottom
    facet!(builder,p7,p8)
    facet!(builder,p8,p9)
    facet!(builder,p9,p10)
    facet!(builder,p10,p11)
    facet!(builder,p11,p12)
    facet!(builder,p12,p7)

    facetregion!(builder,3) # stressor region (bottom, sides, top)
    # bottom
    facet!(builder,p7,p13)
    facet!(builder,p13,p14)
    facet!(builder,p14,p10)

    cellregion!(builder,1) # material 1
    maxvolume!(builder,A_core/6*vol_factor_core)
    regionpoint!(builder,(0,0))

    cellregion!(builder,2) # material 2
    maxvolume!(builder,A_shell/6*vol_factor_shell)
    regionpoint!(builder,(scale[1]+scale[2]/2,0))

    cellregion!(builder,3) # material 3
    maxvolume!(builder,A_stressor/6*vol_factor_stressor)
    regionpoint!(builder,(0,-sqrt(3)/2*d2-δ/2))

    xgrid = simplexgrid(builder)

    hz = 10 * hz_factor
    xgrid = simplexgrid(xgrid,0:hz:scale[4]; bot_offset = 3, top_offset = 6)
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


    ## fix local orderings to avoid negative CellVolumes
    xCellVolumes = xgrid[CellVolumes]
    xCellNodes = xgrid[CellNodes]
    for cell = 1 : num_cells(xgrid)
        if xCellVolumes[cell] < 0
            xCellNodes[[2,1],cell] = xCellNodes[[1,2],cell]
            xCellVolumes[cell] *= -1
        end
    end
    return xgrid

    return xgrid
end

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

    facetregion!(builder,10) # interior boundary
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


    ## fix local orderings to avoid negative CellVolumes
    xCellVolumes = xgrid[CellVolumes]
    xCellNodes = xgrid[CellNodes]
    for cell = 1 : num_cells(xgrid)
        if xCellVolumes[cell] < 0
            xCellNodes[[2,1],cell] = xCellNodes[[1,2],cell]
            xCellVolumes[cell] *= -1
        end
    end
    return xgrid


    return xgrid

end