abstract type MaterialType end
abstract type TestMaterial{x} <: MaterialType where {x <: Float64} end
abstract type GaAs <: MaterialType end
abstract type AlInAs{x} <: MaterialType where {x <: Float64} end
abstract type AlGaAs{x} <: MaterialType where {x <: Float64} end
Base.String(T::Type{<:TestMaterial}) = "Test($(T.parameters[1])))"
Base.String(::Type{GaAs}) = "GaAs"
Base.String(T::Type{<:AlInAs}) = "Ai_$(T.parameters[1])In_$(1 - T.parameters[1])As"
Base.String(T::Type{<:AlGaAs}) = "Ai_$(T.parameters[1])Ga_$(1 - T.parameters[1])As"

struct MaterialDataset{T, MT <: MaterialType}
    ElasticConstants::Dict
    PiezoElectricConstants::Dict
    LatticeConstants::Array{T,1}
    Psp::Array{T,1}
    kappar::T
end

# constructor for isotropic tensor from material data
function IsotropicElasticityTensor(MD::MaterialDataset{T, MT}) where {T,MT}
    # compute Lamé constants in N/(m)^2
    if haskey(d, "λ") && haskey(d, "μ")
        λ = MD("λ")
        μ = MD("μ")
    elseif haskey(d, "C11") && haskey(d, "C44")
        λ = MD("C11") - 2*MD("C44")
        μ = MD("C44")
    else
        @error "material data does not contain enough values to init isotropic elasticity tensor (needs λ, μ or C11, C44)"
    end
    return IsotropicElasticityTensor{3,T}(λ, μ)
end

get_materialtype(::MaterialDataset{T,MT}) where {T,MT} = MT


function MaterialDataset(MT::Type{TestMaterial{x}}) where {x}
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

    return MaterialDataset{Float64, MT}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialDataset(MT::Type{GaAs})

    # elastic constants
    # Phys. Rev. B 95, 245309 (2017)
    # wurtzite
    ElasticConstants = Dict()
    ElasticConstants["C11"] = 1188 # GPa
    ElasticConstants["C12"] =  538 # GPa
    ElasticConstants["C44"] =  594 # GPa

    # piezoelectric constants
    # zinc blende
    PiezoElectricConstants = Dict()
    PiezoElectricConstants["E31wz"] =  0.1328
    PiezoElectricConstants["E33wz"] = -0.2656
    PiezoElectricConstants["E15wz"] = -0.2656 / 6.4825 * 3.9697       # Michele Catti  ASCS2006, Spokane
    PiezoElectricConstants["E14zb"] = PiezoElectricConstants["E15wz"]  # Michele Catti  ASCS2006, Spokane

    # lattice constants
    lc = 5.65325 * ones(Float64,3)
    #lc = [3.9676, 5.65325, 6.4825]

    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant
    kappar = 13.18

    return MaterialDataset{Float64, MT}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialDataset(MT::Type{AlInAs{x}}) where {x}
    # elastic constants
    ElasticConstants = Dict()
    ElasticConstants["C11"] = x*1250 + (1-x)*832.9   # GPa
    ElasticConstants["C12"] = x*534  + (1-x)*452.6   # GPa
    ElasticConstants["C44"] = x*542  + (1-x)*395.9   # GPa

    # piezoelectric constants
    PiezoElectricConstants = Dict()
    PiezoElectricConstants["E31wz"] = x*(0.1)    + (1-x)*0.1            # "Spontaneous polarization and piezoelectric constants of III-V nitrides", Bernadini, Fiorentini and Vanderbilt, Phys Rev B
    PiezoElectricConstants["E33wz"] = x*(-0.01)  + (1-x)*(-0.03)        # "Spontaneous polarization and piezoelectric constants of III-V nitrides", Bernadini, Fiorentini and Vanderbilt, Phys Rev B
    PiezoElectricConstants["E15wz"] = PiezoElectricConstants["E31wz"]   # ?? -- this value for AlInAs not found yet
    PiezoElectricConstants["E14zb"] = x*(-0.048) + (1-x)*(-0.115)       # Phys. Rev. B 84, 195207 (2011)

    # lattice constants
    lc = (x*5.6605 + (1-x)*6.0553) * ones(Float64,3)

    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant
    kappar = 12.605

    return MaterialDataset{Float64, MT}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end

function MaterialDataset(MT::Type{AlGaAs{x}}) where {x}
    # elastic constants
    ElasticConstants = Dict()
    ElasticConstants["C11"] = x*1250 + (1-x)*1188    # GPa
    ElasticConstants["C12"] = x*534  + (1-x)*538     # GPa
    ElasticConstants["C44"] = x*542  + (1-x)*594     # GPa

    # piezoelectric constants (take from GaAs / AlAs data above)
    PiezoElectricConstants = Dict()
    PiezoElectricConstants["E31wz"] = x*(0.1)    + (1-x)*(0.1328)
    PiezoElectricConstants["E33wz"] = x*(-0.01)  + (1-x)*(-0.2656)
    PiezoElectricConstants["E15wz"] = x*(0.1)    + (1-x)*(-0.2656 / 6.4825 * 3.9697)
    PiezoElectricConstants["E14zb"] = x*(-0.048) + (1-x)*(-0.2656 / 6.4825 * 3.9697)

    # lattice constants
    lc = (x*5.6605 + (1-x)*5.65325) * ones(Float64,3)
    
    # C/m2 spontaneous polarization
    Psp = zeros(Float64,3)

    # relative dielectric constant 
    kappar = 12.605 # NEED TO CORRECT FROM LITERATURE
    
    return MaterialDataset{Float64, MT}(ElasticConstants, PiezoElectricConstants, lc, Psp, kappar)
end


struct MaterialData{T, MST <: MaterialStructureType}
    data::Array{MaterialDataset,1}
    TensorC::Array{<:ElasticityTensorType{T},1}         # elastitity tensor C
    TensorE::Array{<:PiezoElectricityTensorType{T},1}   # piezo-electricity tensor E
end

### coefficients ###
function set_data(materials::Array{DataType,1}, MST::Type{<:MaterialStructureType})

    data = Array{MaterialDataset,1}(undef, length(materials))
    C = Array{CustomMatrixElasticityTensor{Float64},1}(undef,length(materials))
    E = Array{CustomMatrixPiezoElectricityTensor{Float64},1}(undef,length(materials))

    for m = 1 : length(materials)

        ## get material data set
        data[m] = MaterialDataset(materials[m])

        ## setup tensors
        if MST <: ZincBlende001
            # stiffness tensor in N/(nm)^2 # todo: move scaling to solver
            # Equation (4)
            # in "Symmetry-adapted calculations of strain and polarization fields in (111)-oriented 
            # zinc-blende quantum dots" by Schulz et. al
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            C[m] = CustomMatrixElasticityTensor(
              1e-9*[ C11 C12 C12   0   0   0
                     C12 C11 C12   0   0   0
                     C12 C12 C11   0   0   0
                       0   0   0 C44   0   0
                       0   0   0   0 C44   0
                       0   0   0   0   0 C44 ])
            # piezoelectric tensor in C/(m^2)
            # Equation (22)
            E14zb = data[m].PiezoElectricConstants["E14zb"]
            E[m] = CustomMatrixPiezoElectricityTensor(
                  [0 0 0 E14zb     0     0
                   0 0 0     0 E14zb     0
                   0 0 0     0     0 E14zb])
        elseif MST <: ZincBlende2D
            # stiffness tensor in N/(nm)^2
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            C[m] = CustomMatrixElasticityTensor(
                   1e-9*[ C11 C12   0
                          C12 C11   0
                            0   0 C44 ])
            # piezoelectric tensor in C/(m^2)
            E14zb = data[m].PiezoElectricConstants["E14zb"]
            E[m] = CustomMatrixPiezoElectricityTensor(
                [E14zb     0 0 0 0 0
                     0 E14zb 0 0 0 0])
        elseif MST <: ZincBlende111
            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            sr2 = sqrt(2)
            C11p = (1/2)*(C11 + C12) + C44
            C12p = (1/6)*(C11 + 5*C12) - (1/3) * C44
            C44p = (1/3)*(C11 + C12) + (1/3) * C44
            C33p = (3/2)*C11p - (1/2)*C12p -C44p
            C13p = -(1/2)*C11p + (3/2)*C12p + C44p
            C15p = (1/sr2)*C11p - (1/sr2)*C12p - sr2*C44p
            C66p = (1/2)*(C11p - C12p)

            # stiffness tensor in N/(nm)^2
            # Equation (14)
            C[m] = CustomMatrixElasticityTensor(
              1e-9*[ C11p  C12p C13p     0  C15p     0
                     C12p  C11p C13p     0 -C15p     0
                     C13p  C13p C33p     0     0     0
                        0     0    0  C44p     0 -C15p
                     C15p -C15p    0     0  C44p     0
                        0     0    0 -C15p     0  C66p])

            # piezoelectric tensor in C/(m^2)
            # Equation (27)

            E14zb = data[m].PiezoElectricConstants["E14zb"]
            E11 = - sqrt(2/3) * E14zb
            E12 = sqrt(2/3) * E14zb
            E15 = - sqrt(1/3) * E14zb
            E31 = E15
            E33 = 2/sqrt(3) * E14zb

            E[m] = CustomMatrixPiezoElectricityTensor(
                   [E11 E12 0   0   E15 0
                    0   0   0   E15 0   E12
                    E31 E31 E33 0   0   0])
        elseif MST <: Wurtzite0001

            # stiffness tensor in N/(nm)^2
            # Equation (7)

            C11 = data[m].ElasticConstants["C11"]
            C12 = data[m].ElasticConstants["C12"]
            C44 = data[m].ElasticConstants["C44"]
            C11wz = (1/6) * (3*C11 + 3 * C12 + 6 * C44)
            C33wz = (1/6) * (2*C11 + 4 * C12 + 8 * C44)
            C12wz = (1/6) * (1*C11 + 5 * C12 - 2 * C44)
            C13wz = (1/6) * (2*C11 + 4 * C12 - 4 * C44)
            C44wz = (1/6) * (2*C11 - 2 * C12 + 2 * C44)
            C66wz = (1/6) * (1*C11 - 1 * C12 + 4 * C44)
            C[m] = CustomMatrixElasticityTensor(
                   1e-9*[ C11wz C12wz C13wz 0     0     0
                          C12wz C11wz C13wz 0     0     0
                          C13wz C13wz C33wz 0     0     0
                          0     0     0     C44wz 0     0
                          0     0     0     0     C44wz 0
                          0     0     0     0     0     C66wz ])


            # piezoelectric tensor in C/(m^2)
            # Equation (25)
            E31wz = data[m].PiezoElectricConstants["E31wz"]
            E33wz = data[m].PiezoElectricConstants["E33wz"]
            E15wz = data[m].PiezoElectricConstants["E15wz"]
            E[m] = CustomMatrixPiezoElectricityTensor(
                   [     0     0     0     0  E15wz  0 
                         0     0     0 E15wz      0  0 
                     E31wz E31wz E33wz     0      0  0 ])
        end
    end

    return MaterialData{Float64,MST}(data,C,E)
end
