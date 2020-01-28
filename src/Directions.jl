"""
# module Directions



# Examples

```jldoctest
julia>
```
"""
module Directions

using StaticArrays: FieldVector

using Crystallography: AbstractSpace, RealSpace, ReciprocalSpace

export MillerIndices, MillerBravaisIndices

struct MillerIndices{S<:AbstractSpace,T<:Integer} <: FieldVector{3,T}
    i::T
    j::T
    k::T
    function MillerIndices{S,T}(i, j, k) where {S,T}
        x = [i, j, k]
        i, j, k = iszero(x) ? x : x .÷ gcd(x)
        return new(i, j, k)
    end
end
MillerIndices{S}(i::T, j::T, k::T) where {S,T} = MillerIndices{S,T}(i, j, k)
MillerIndices{S}(x::AbstractVector) where {S} = MillerIndices{S}(x...)
MillerIndices{S}(x::Tuple) where {S} = MillerIndices{S}(collect(x))

struct MillerBravaisIndices{S<:AbstractSpace,T<:Integer} <: FieldVector{4,T}
    i::T
    j::T
    k::T
    l::T
    function MillerBravaisIndices{S,T}(i, j, k, l) where {S,T}
        x = [i, j, k, l]
        i, j, k, l = iszero(x) ? x : x .÷ gcd(x)
        return new(i, j, k, l)
    end
end
MillerBravaisIndices{S}(i::T, j::T, k::T, l::T) where {S,T} =
    MillerBravaisIndices{S,T}(i, j, k, l)
MillerBravaisIndices{S}(x::AbstractVector) where {S} = MillerBravaisIndices{S}(x...)
MillerBravaisIndices{S}(x::Tuple) where {S} = MillerBravaisIndices{S}(collect(x))

Base.getproperty(x::MillerIndices{RealSpace}, name::Symbol) =
    getfield(x, Dict(:u => :i, :v => :j, :w => :k)[name])
Base.getproperty(x::MillerIndices{ReciprocalSpace}, name::Symbol) =
    getfield(x, Dict(:h => :i, :k => :j, :l => :k)[name])
Base.getproperty(x::MillerBravaisIndices{RealSpace}, name::Symbol) =
    getfield(x, Dict(:u => :i, :v => :j, :t => :k, :w => :l)[name])
Base.getproperty(x::MillerBravaisIndices{ReciprocalSpace}, name::Symbol) =
    getfield(x, Dict(:h => :i, :k => :j, :i => :k, :l => :l)[name])

Base.convert(::Type{<:MillerIndices}, mb::MillerBravaisIndices) =
    error("Unsupported operation!")
Base.convert(::Type{<:MillerBravaisIndices}, m::MillerIndices) =
    error("Unsupported operation!")
Base.convert(T::Type{<:MillerIndices}, m::MillerIndices) =
    isa(m, T) ? m : error("Unsupported operation!")
Base.convert(T::Type{<:MillerBravaisIndices}, mb::MillerBravaisIndices) =
    isa(mb, T) ? mb : error("Unsupported operation!")
Base.convert(::Type{MillerIndices{T}}, mb::MillerBravaisIndices{T}) where {T<:RealSpace} =
    MillerIndices{T}(2 * mb.u + mb.v, 2 * mb.v + mb.u, mb.w)
Base.convert(
    ::Type{MillerIndices{T}},
    mb::MillerBravaisIndices{T},
) where {T<:ReciprocalSpace} = MillerIndices{T}(mb.h, mb.k, mb.l)
Base.convert(::Type{MillerBravaisIndices{T}}, m::MillerIndices{T}) where {T<:RealSpace} =
    MillerBravaisIndices{T}(2 * m.u - m.v, 2 * m.v - m.u, -(m.u + m.v), 3 * m.w)
Base.convert(
    ::Type{MillerBravaisIndices{T}},
    m::MillerIndices{T},
) where {T<:ReciprocalSpace} = MillerBravaisIndices{T}(m.h, m.k, -(m.h + m.k), m.l)

end
