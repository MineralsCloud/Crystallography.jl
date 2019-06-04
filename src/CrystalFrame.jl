"""
# module CrystalFrame



# Examples

```jldoctest
julia>
```
"""
module CrystalFrame

using Crystallography

export AbstractCoordinatesTransformation,
    CrystalFromCartesian,
    CartesianFromCrystal,
    transformation

abstract type AbstractCoordinatesTransformation{T <: BravaisLattice} end
struct CrystalFromCartesian{T} <: AbstractCoordinatesTransformation{T} end
struct CartesianFromCrystal{T} <: AbstractCoordinatesTransformation{T} end

transformation(::Type{T}, a::CrystalCoordinates) where {T <: BravaisLattice} = transformation(CartesianFromCrystal{T}, collect(a))
transformation(::Type{T}, a::CartesianCoordinates) where {T <: BravaisLattice} = transformation(CrystalFromCartesian{T}, collect(a))

end
