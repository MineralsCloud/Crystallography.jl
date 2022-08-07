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
LatticeSystem
Bravais
Lattice
basis_vectors
latticesystem
latticeconstants
periodicity
supercell
```

### Reciprocal space

Note that we take ``2\pi`` as ``1``, not the solid-state physics convention.

```@docs
ReciprocalPoint
ReciprocalLattice
MonkhorstPackGrid
reciprocal
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
crystaldensity
```
