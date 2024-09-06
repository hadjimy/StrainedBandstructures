# Material Data

## MaterialParameters

This structure stores all implemented parameters for a material.

```@docs
MaterialParameters
```

## Material structure types

```@autodocs
Modules = [StrainedBandstructures]
Pages = ["materialstructuretype.jl"]
Order   = [:type, :function]
```

## Available Materials

```@docs
MaterialType
GaAs
AlInAs
AlGaAs
InGaAs
```

## HeteroStructureData

A hetero-structure consists of several materials related to the region numbers of the
grid. The HeteroStructureData is a structure that stores the MaterialParameters
of all materials/regions and the derived elasticity and piezo-eletricity tensors.

```@autodocs
Modules = [StrainedBandstructures]
Pages = ["heterostructures.jl"]
Order   = [:type, :function]
```

