module NanoWiresJulia

using GradientRobustMultiPhysics
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using SimplexGridFactory
using Triangulate
using TetGen
using LinearAlgebra
using Printf



## problem: loading and saving grids leads to ElementGeometries -> DataType conversion (by DrWarson/JLD2?) which has to be reverted
## after loading (until this is fixed ina more elegant way)
function repair_grid!(xgrid::ExtendableGrid)
    xgrid.components[CellGeometries] = VectorOfConstants{ElementGeometries,Int}(xgrid.components[CellGeometries][1], num_cells(xgrid)) 
    xgrid.components[FaceGeometries] = VectorOfConstants{ElementGeometries,Int}(xgrid.components[FaceGeometries][1], length(xgrid.components[FaceGeometries])) 
    xgrid.components[BFaceGeometries] = VectorOfConstants{ElementGeometries,Int}(xgrid.components[BFaceGeometries][1], length(xgrid.components[BFaceGeometries])) 

    xgrid.components[UniqueCellGeometries] = Vector{ElementGeometries}([xgrid.components[CellGeometries][1]])
    xgrid.components[UniqueFaceGeometries] = Vector{ElementGeometries}([xgrid.components[FaceGeometries][1]])
    xgrid.components[UniqueBFaceGeometries] = Vector{ElementGeometries}([xgrid.components[BFaceGeometries][1]])
end
export repair_grid!

include("grids.jl")
export bimetal_strip3D, bimetal_strip2D, bimetal_strip3D_middle_layer, bimetal_tensorgrid, bimetal_tensorgrid_uniform
export condensator3D, condensator3D_tensorgrid, condensator2D
export nonpolarquantumwell3D, nonpolarquantumwell2D
export nanowire_grid, nanowire_tensorgrid

include("materialstructuretype.jl")
export MaterialStructureType
export ZincBlende2D, ZincBlende001, ZincBlende111, Wurtzite0001

include("tensors.jl")
export ElasticityTensorType, PiezoElectricityTensorType
export IsotropicElasticityTensor
export CustomMatrixElasticityTensor, CustomMatrixPiezoElectricityTensor
export apply_elasticity_tensor!, apply_piezoelectricity_tensor!

include("materials.jl")
export MaterialType
export TestMaterial, GaAs, AlInAs
export MaterialDataset
export get_materialtype
export MaterialData
export set_data

include("strains.jl")
export StrainType
export LinearStrain, LinearStrain2D, LinearStrain3D
export NonlinearStrain, NonlinearStrain2D, NonlinearStrain3D
export eval_strain!

include("pdeoperators.jl")
export get_displacement_operator, get_displacement_operator_new
export get_polarisation_from_strain_operator
export get_polarisation_laplacian_operator
export get_polarisation_operator # wip
export get_energy_integrator

include("plane_cuts.jl")
export tet_x_plane!, plane_cut, get_cutgrids
export perform_simple_plane_cuts, perform_plane_cuts

include("export.jl")
export exportVTK

include("statistics.jl")
export compute_statistics

include("solvers.jl")
export solve_by_embedding
export solve_by_embedding!
export solve_by_damping
export solve_lowlevel

greet() = print("Hello World!")

end # module
