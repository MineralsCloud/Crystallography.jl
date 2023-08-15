# Public API

## Contents

```@contents
Pages = ["public.md"]
Depth = 2
```

## Index

```@index
Pages = ["public.md"]
```

## Public interface

### Lattice

```@docs
CrystalSystem
Triclinic
Monoclinic
Orthorhombic
Tetragonal
Cubic
Trigonal
Hexagonal
Centering
BaseCentering
Primitive
BodyCentering
FaceCentering
RhombohedralCentering
BaseCentering
Bravais
Lattice
centering
crystalsystem
basis_vectors
cellparameters
supercell
```

### Reciprocal space

Note that we take ``2pi`` as ``1``, not the solid-state physics convention.

```@docs
ReciprocalPoint
ReciprocalLattice
inv
reciprocal_mesh
coordinates
weights
```

### Miller and Miller–Bravais indices

```@docs
Miller
MillerBravais
ReciprocalMiller
ReciprocalMillerBravais
family
@m_str
```

### Metric tensor

```@docs
MetricTensor
directioncosine
directionangle
distance
interplanar_spacing
```

### Transformations

```@docs
CartesianFromFractional
FractionalFromCartesian
PrimitiveFromStandardized
StandardizedFromPrimitive
```

### Others

```@docs
cellvolume
```

### Symmetry

```@docs
SeitzOperator
```
