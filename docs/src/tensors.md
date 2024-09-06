# Tensors

## Custom tensors 

Custom tensors are based on a matrix representation that maps one vector to another, e.g. a Voigt representation of the strain to the Voigt representation of the stress. The matrix representation can be defined as a DenseMatrix or a SparseMatrix.

```@docs
AbstractTensor
SparseMatrixTensor
DenseMatrixTensor
```

Currently DenseMatrixTensors are used to represent the elasticity and piezoeletricity tensors.

## Isotropic Elasticity tensors

For linear elasticity, a special constructor for the isotropic elasticity tensor is available:

```@docs
IsotropicElasticityTensor
IsotropicElasticityTensor(MD::MaterialParameters{T, MT, MST}) where {T,MT,MST}
```