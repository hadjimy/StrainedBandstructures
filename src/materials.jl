
"""
    $(TYPEDEF)
"""
abstract type MaterialType end
abstract type TestMaterial{x} <: MaterialType where {x <: Float64} end
"""
    $(TYPEDEF)

    Gallium Arsenite
"""
abstract type GaAs <: MaterialType end
"""
    $(TYPEDEF)

    Alluminum-Indium-Arsenite compound where x ...
"""
abstract type AlInAs{x} <: MaterialType where {x <: Float64} end
"""
    $(TYPEDEF)

    Alluminum-Gallium-Arsenite compound where x ...
"""
abstract type AlGaAs{x} <: MaterialType where {x <: Float64} end
"""
    $(TYPEDEF)

    Indium-Gallium-Arsenite compound where x ...
"""
abstract type InGaAs{x} <: MaterialType where {x <: Float64} end

Base.String(T::Type{<:TestMaterial}) = "Test($(T.parameters[1])))"
Base.String(::Type{GaAs}) = "GaAs"
Base.String(T::Type{<:AlInAs}) = "Al_$(T.parameters[1])In_$(1 - T.parameters[1])As"
Base.String(T::Type{<:AlGaAs}) = "Al_$(T.parameters[1])Ga_$(1 - T.parameters[1])As"
Base.String(T::Type{<:InGaAs}) = "In_$(T.parameters[1])Ga_$(1 - T.parameters[1])As"

"""
    $(TYPEDEF)

Data set that collects all parameters for a material (region), i.e.
elastic constants, piezoelextric constants, lattice constants, the spontaneous polarization constant and the relative dielectric constant.
"""
struct MaterialParameters{T, MT <: MaterialType, MST <: MaterialStructureType}
    """
    Elastic Constants
    """
    ElasticConstants::Dict
    """
    Piezoelectric constants
    """
    PiezoElectricConstants::Dict
    """
    Lattice constants
    """
    LatticeConstants::Array{T,1}
    """
    spontaneous polarization
    """
    Psp::Array{T,1}
    """
    relative dielectric constant 
    """
    kappar::T
end


function get_kappar(mds::MaterialParameters)
    MT=get_materialtype(mds)
    return MT(mds.kappar)
end

get_materialtype(::MaterialParameters{T,MT,MST}) where {T,MT,MST} = MT


"""
````
function MaterialParameters(MT::MaterialType)
````

Returns the MaterialDateset for MaterialType MT.

"""
function MaterialParameters(MT::Type{TestMaterial{x}}) where {x}
    # elastic constants
    ElasticConstants = Dict()
    ElasticConstants["C11"] = 1200 # GPa
    ElasticConstants["C12"] =  300 # GPa
    ElasticConstants["C44"] =  500 # GPa

    # piezoelectric constants
    PiezoElectricConstants = Dict()
    PiezoElectricConstants["E31wz"] =  0.05
    PiezoElectricConstants["E33wz"] = -0.02
    PiezoElectricConstants["E14zb"] = -0.08

    # lattice constants
    lc = 5*(1+x) * ones(Float64,3)
    
    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant 
    kappar = 10.0

    return MaterialParameters{Float64, MT}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialParameters(MT::Type{GaAs}, MST::Type{<:MaterialStructureType})

    # elastic constants
    # Phys. Rev. B 95, 245309 (2017)
    # wurtzite
    ElasticConstants = Dict() # from Vurgaftman, see also SPHInX GaAs.sx file
    ElasticConstants["C11"] = 122.1 # GPa (118.8)
    ElasticConstants["C12"] =  56.6 # GPa (53.8)
    ElasticConstants["C44"] =  60.0 # GPa (59.4)

    # piezoelectric constants
    # zinc blende
    PiezoElectricConstants = Dict()
    # Quasi-cubic approximation using the values from Beya-Wakata et al., Phys. Rev. B 84, 195207 (2011)
    # see also Bernardini et al., Phys. Rev. B 56, R10024 (1997)
    PiezoElectricConstants["E14zb"] = -0.2381 # Beya-Wakata et al., Phys. Rev. B 84, 195207 (2011), other value from O.M.: PiezoElectricConstants["E15wz"] # Michele Catti  ASCS2006, Spokane
    PiezoElectricConstants["E31wz"] = -1/sqrt(3)*PiezoElectricConstants["E14zb"] # other value from O.M.: 0.1328
    PiezoElectricConstants["E33wz"] = 2/sqrt(3)*PiezoElectricConstants["E14zb"]  # other value from O.M.: -0.2656
    PiezoElectricConstants["E15wz"] = PiezoElectricConstants["E31wz"] # other value from O.M.: -0.2656 / 6.4825 * 3.9697 # Michele Catti  ASCS2006, Spokane

    # lattice constants
    if MST <: ZincBlende001
        lc = 5.6532 * ones(Float64,3) # 5.750 from https://materialsproject.org/materials/mp-2534?formula=GaAs (space group: F4¯3m)
    elseif MST <: ZincBlende111_C14 || MST <: ZincBlende111_C15 || MST <: ZincBlende111_C14_C15
        lc = [5.6532/sqrt(2), 5.6532/sqrt(2), 5.6532/sqrt(3/4)] # From GaAs.sx file in SPHInX
    elseif  MST <: Wurtzite0001
        lc = [3.994, 3.994, 6.586] # from https://materialsproject.org/materials/mp-8883?formula=GaAs (space group: P6₃mc)
    end

    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant
    kappar = 12.9 # http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html

    return MaterialParameters{Float64, MT, MST}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialParameters(MT::Type{AlInAs{x}}, MST::Type{<:MaterialStructureType}) where {x}
    # elastic constants
    ElasticConstants = Dict() # from Vurgaftman, see also SPHInX AlInAs.sx file
    ElasticConstants["C11"] = x*125.0 + (1-x)*83.29 # GPa
    ElasticConstants["C12"] = x*53.4  + (1-x)*45.26 # GPa
    ElasticConstants["C44"] = x*54.2  + (1-x)*39.59 # GPa

    # piezoelectric constants
    PiezoElectricConstants = Dict()
    # Quasi-cubic approximation using the values from Beya-Wakata et al., Phys. Rev. B 84, 195207 (2011)
    # see also Bernardini et al., Phys. Rev. B 56, R10024 (1997)
    PiezoElectricConstants["E14zb"] = x*(-0.048) + (1-x)*(-0.115)                # Phys. Rev. B 84, 195207 (2011)
    PiezoElectricConstants["E31wz"] = -1/sqrt(3)*PiezoElectricConstants["E14zb"] # "Spontaneous polarization and piezoelectric constants of III-V nitrides", Bernadini, Fiorentini and Vanderbilt, Phys Rev B
    PiezoElectricConstants["E33wz"] = 2/sqrt(3)*PiezoElectricConstants["E14zb"]  # "Spontaneous polarization and piezoelectric constants of III-V nitrides", Bernadini, Fiorentini and Vanderbilt, Phys Rev B
    PiezoElectricConstants["E15wz"] = PiezoElectricConstants["E31wz"]

    # lattice constants
    if MST <: ZincBlende001
        lc = (x*5.6611 + (1-x)*6.0583) * ones(Float64,3) # x*5.676 + (1-x)*6.107 from https://materialsproject.org/materials/mp-2172?formula=AlAs and https://materialsproject.org/materials/mp-20305?formula=InAs (space groups: F4¯3m)
    elseif MST <: ZincBlende111_C14 || MST <: ZincBlende111_C15 || MST <: ZincBlende111_C14_C15
        lc = [(x*5.6611 + (1-x)*6.0583)/sqrt(2), (x*5.6611 + (1-x)*6.0583)/sqrt(2), (x*5.6611 + (1-x)*6.0583)/sqrt(3/4)] # From AlInAs.sx file in SPHInX
    elseif  MST <: Wurtzite0001
        lc = [(x*4.002 + (1-x)*4.311), (x*4.002 + (1-x)*4.311), (x*6.582 + (1-x)*7.093)] # from https://materialsproject.org/materials/mp-8881?formula=AlAs and https://materialsproject.org/materials/mp-1007652?formula=InAs (space group: P6₃mc)
    end

    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant
    kappar = x*10.06  + (1-x)*15.15 # from https://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html and https://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html

    return MaterialParameters{Float64, MT, MST}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialParameters(MT::Type{AlGaAs{x}}, MST::Type{<:MaterialStructureType}) where {x}
    # elastic constants (taken from GaAs / AlAs data above)
    ElasticConstants = Dict()
    ElasticConstants["C11"] = x*125.0 + (1-x)*122.1 # GPa
    ElasticConstants["C12"] = x*53.4  + (1-x)*56.6  # GPa
    ElasticConstants["C44"] = x*54.2  + (1-x)*60.0  # GPa

    # piezoelectric constants (taken from GaAs / AlAs data above)
    PiezoElectricConstants = Dict()
    PiezoElectricConstants["E31wz"] = x*(0.1)    + (1-x)*(0.1328)
    PiezoElectricConstants["E33wz"] = x*(-0.01)  + (1-x)*(-0.2656)
    PiezoElectricConstants["E15wz"] = x*(0.1)    + (1-x)*(-0.2656 / 6.4825 * 3.9697)
    PiezoElectricConstants["E14zb"] = x*(-0.048) + (1-x)*(-0.2656 / 6.4825 * 3.9697)

    # lattice constants (taken from GaAs / AlAs data above)
    if MST <: ZincBlende001
        lc = (x*5.6611 + (1-x)*5.6532) * ones(Float64,3)
    elseif MST <: ZincBlende111_C14 || MST <: ZincBlende111_C15 || MST <: ZincBlende111_C14_C15
        lc = [(x*5.6611 + (1-x)*5.6532)/sqrt(2), (x*5.6611 + (1-x)*5.6532)/sqrt(2), (x*5.6611 + (1-x)*5.6532)/sqrt(3/4)]
    elseif  MST <: Wurtzite0001
        lc = [(x*4.002 + (1-x)*3.994), (x*4.002 + (1-x)*3.994), (x*6.582 + (1-x)*6.586)]
    end
    
    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant (taken from GaAs / AlAs data above)
    kappar = x*10.06  + (1-x)*12.9

    return MaterialParameters{Float64, MT, MST}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialParameters(MT::Type{InGaAs{x}}, MST::Type{<:MaterialStructureType}) where {x}
    # elastic constants (taken from InAs / GaAs data above)
    ElasticConstants = Dict()
    ElasticConstants["C11"] = x*83.29 + (1-x)*122.1 # GPa
    ElasticConstants["C12"] = x*45.26 + (1-x)*56.6  # GPa
    ElasticConstants["C44"] = x*39.59 + (1-x)*60.0  # GPa

    # piezoelectric constants (taken from InAs / GaAs data above)
    PiezoElectricConstants = Dict()
    PiezoElectricConstants["E31wz"] = x*(0.1)    + (1-x)*(0.1328)
    PiezoElectricConstants["E33wz"] = x*(-0.03)  + (1-x)*(-0.2656)
    PiezoElectricConstants["E15wz"] = x*(0.1)    + (1-x)*(-0.2656 / 6.4825 * 3.9697)
    PiezoElectricConstants["E14zb"] = x*(-0.115) + (1-x)*(-0.2656 / 6.4825 * 3.9697)

    # lattice constants (taken from InAs / GaAs data above)
    if MST <: ZincBlende001
        lc = (x*6.0583 + (1-x)*5.6532) * ones(Float64,3)
    elseif MST <: ZincBlende111_C14 || MST <: ZincBlende111_C15 || MST <: ZincBlende111_C14_C15
        lc = [(x*6.0583 + (1-x)*5.6532)/sqrt(2), (x*6.0583 + (1-x)*5.6532)/sqrt(2), (x*6.0583 + (1-x)*5.6532)/sqrt(3/4)]
    elseif  MST <: Wurtzite0001
        lc = [(x*4.311 + (1-x)*3.994), (x*4.311 + (1-x)*3.994), (x*7.093 + (1-x)*6.586)]
    end

    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant (taken from InAs / GaAs data above)
    kappar = x*15.15  + (1-x)*12.9

    return MaterialParameters{Float64, MT, MST}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

