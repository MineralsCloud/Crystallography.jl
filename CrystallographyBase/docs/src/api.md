```@meta
CurrentModule = CrystallographyBase
```

# API

```@contents
Pages = ["api.md"]
Depth = 3
```

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

```@docs
ReciprocalPoint
ReciprocalLattice
inv
reciprocal_mesh
coordinates
weights
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
