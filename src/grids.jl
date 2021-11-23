function bimetal_strip3D(; material_border = 0.5, scale = [1,1,1], anisotropy = [1,1,2], reflevel = 1, maxvol = prod(scale./anisotropy)/8)

    scale ./= anisotropy

    builder=SimplexGridBuilder(Generator=TetGen)

    @info "Generating bimetal grid for scale = $scale"

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

function bimetal_strip2D(; material_border = 0.5, scale = [1,1], anisotropy = [1,1], reflevel = 1, maxvol = prod(scale./anisotropy)/4)

    scale ./= anisotropy

    builder=SimplexGridBuilder(Generator=Triangulate)

    @info "Generating 2d bimetal grid for scale = $scale"

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

    p1=point!(builder,0,0,0)                                                
    p2=point!(builder,B,0,0)                                         
    p3=point!(builder,B,0,A)                                 
    p4=point!(builder,B,0,A+d)    
    p5=point!(builder,B,0,2*A+d)                                                
    p6=point!(builder,0,0,2*A+d)              
    p7=point!(builder,0,0,A+d)                           
    p8=point!(builder,0,0,A)     

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
    facetregion!(builder,2) # left and right boundary of bottom plate
    facet!(builder,p2,p3)
    facet!(builder,p8,p1)

    facetregion!(builder,2) # gap boundary left and right
    facet!(builder,p3 ,p4)
    facet!(builder,p7 ,p8)

    facetregion!(builder,3) # boundary of top plate (left, top, right)
    facet!(builder,p4, p5)
    facet!(builder,p5, p6)
    facet!(builder,p6, p7)

    facetregion!(builder,99) # interior facets to split materials
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

