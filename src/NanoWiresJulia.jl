module NanoWiresJulia

using GradientRobustMultiPhysics
using ExtendableGrids
using GridVisualize
using SimplexGridFactory
using Triangulate
using TetGen
using LinearAlgebra
using Printf

include("grids.jl")
export bimetal_strip3D, bimetal_strip2D
export condensator3D, condensator2D

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
export get_displacement_operator
export get_polarisation_from_strain_operator
export get_polarisation_laplacian_operator
export get_polarisation_operator # wip
export get_energy_integrator

include("plane_cuts.jl")
export tet_x_plane!, plane_cut, get_cutgrids, perform_plane_cuts

include("export.jl")
export writeVTK

include("statistics.jl")
export compute_statistics

include("solvers.jl")
export solve_by_embedding
export solve_by_damping

greet() = print("Hello World!")

end # module
