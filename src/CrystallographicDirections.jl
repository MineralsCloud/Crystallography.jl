"""
# module CrystallographicDirections



# Examples

```jldoctest
julia>
```
"""
module CrystallographicDirections

using StaticArrays: FieldVector

using Crystallography

export MillerIndices,
    MillerBravaisIndices

struct MillerIndices{S,T <: Integer} <: FieldVector{3,T}
    i::T
    j::T
    k::T
    function MillerIndices{S,T}(i, j, k) where {S <: AbstractSpace,T <: Integer}
        x = [i, j, k]
        i, j, k = iszero(x) ? x : x .÷ gcd(x)
        new(i, j, k)
    end
end
MillerIndices{S}(i::T, j::T, k::T) where {S <: AbstractSpace,T <: Integer} = MillerIndices{S,T}(i, j, k)
MillerIndices{S}(x::AbstractVector{T}) where {S <: AbstractSpace,T <: Integer} = MillerIndices{S}(x...)
MillerIndices{S}(x::Tuple) where {S} = MillerIndices{S}(collect(x))

struct MillerBravaisIndices{S,T <: Integer} <: FieldVector{4,T}
    i::T
    j::T
    k::T
    l::T
    function MillerBravaisIndices{S,T}(i, j, k, l) where {S <: AbstractSpace,T <: Integer}
        x = [i, j, k, l]
        i, j, k, l = iszero(x) ? x : x .÷ gcd(x)
        new(i, j, k, l)
    end
end
MillerBravaisIndices{S}(i::T, j::T, k::T, l::T) where {S <: AbstractSpace,T <: Integer} = MillerBravaisIndices{S,T}(i, j, k, l)
MillerBravaisIndices{S}(x::AbstractVector{T}) where {S <: AbstractSpace,T <: Integer} = MillerBravaisIndices{S}(x...)
MillerBravaisIndices{S}(x::Tuple) where {S} = MillerBravaisIndices{S}(collect(x))

function Base.getproperty(x::MillerIndices{RealSpace}, name::Symbol)
    getfield(x, Dict(:u => :i, :v => :j, :w => :k)[name])
end
function Base.getproperty(x::MillerIndices{ReciprocalSpace}, name::Symbol)
    getfield(x, Dict(:h => :i, :k => :j, :l => :k)[name])
end
function Base.getproperty(x::MillerBravaisIndices{RealSpace}, name::Symbol)
    getfield(x, Dict(:u => :i, :v => :j, :t => :k, :w => :l)[name])
end
function Base.getproperty(x::MillerBravaisIndices{ReciprocalSpace}, name::Symbol)
    getfield(x, Dict(:h => :i, :k => :j, :i => :k, :l => :l)[name])
end

Base.convert(::Type{<: MillerIndices}, mb::MillerBravaisIndices{RealSpace}) = MillerIndices{RealSpace}(2mb.u + mb.v, 2mb.v + mb.u, mb.w)
Base.convert(::Type{<: MillerIndices}, mb::MillerBravaisIndices{ReciprocalSpace}) = MillerIndices{ReciprocalSpace}(mb.h, mb.k, mb.l)
Base.convert(::Type{<: MillerBravaisIndices}, m::MillerIndices{RealSpace}) = MillerBravaisIndices{RealSpace}(2m.u - m.v, 2m.v - m.u, -(m.u + m.v), 3m.w)
Base.convert(::Type{<: MillerBravaisIndices}, m::MillerIndices{ReciprocalSpace}) = MillerBravaisIndices{ReciprocalSpace}(m.h, m.k, -(m.h + m.k), m.l)

end
