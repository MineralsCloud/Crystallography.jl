# Public API

## Contents

```@contents
Pages = ["public.md"]
Depth = 3
```

## Index

```@index
Pages = ["public.md"]
```

## Public interface

### Lattice and Cell

```@docs
AbstractLattice
Lattice
isrighthanded
basisvectors
latticesystem
latticeconstants
periodicity
Cell
supercell
eachatom
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

### Metric tensor

```@docs
MetricTensor
distance
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
