export Miller, MillerBravais, RealSpace, ReciprocalSpace
export @m_str, @mb_str

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

struct Miller{S<:AbstractSpace} <: AbstractVector{Int}
    data::SVector{3,Int}
    Miller{S}(v) where {S} = new(iszero(v) ? v : v .÷ gcd(v))
end
Miller{S}(i, j, k) where {S} = Miller{S}([i, j, k])

struct MillerBravais{S<:AbstractSpace} <: AbstractVector{Int}
    data::SVector{4,Int}
    function MillerBravais{S}(v) where {S}
        @assert(
            v[3] == -v[1] - v[2],
            "the 3rd index of `MillerBravais` should equal to the negation of the first two!"
        )
        return new(iszero(v) ? v : v .÷ gcd(v))
    end
end
MillerBravais{S}(i, j, k, l) where {S} = MillerBravais{S}([i, j, k, l])

# This is a helper type and should not be exported!
const INDICES = Union{Miller,MillerBravais}

# This is a helper function and should not be exported!
function _indices_str(r::Regex, s::AbstractString, ::Type{T}) where {T<:INDICES}
    m = match(r, strip(s))
    if m === nothing
        error("not a valid expression!")
    else
        brackets = first(m.captures) * last(m.captures)
        x = (parse(Int, x) for x in m.captures[2:(end-1)])
        if brackets ∈ ("()", "{}")
            return T{ReciprocalSpace}(x...)
        elseif brackets ∈ ("[]", "<>")
            return T{RealSpace}(x...)
        else
            error("not a valid expression!")
        end
    end
end

macro m_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    _indices_str(r, s, Miller)
end

macro mb_str(s)
    r =
        r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    _indices_str(r, s, MillerBravais)
end

Base.size(::Miller) = (3,)
Base.size(::MillerBravais) = (4,)

Base.IndexStyle(::Type{<:Union{Miller,MillerBravais}}) = IndexLinear()

Base.getindex(x::Union{Miller,MillerBravais}, i::Integer) = getindex(x.data, i)

Base.convert(::Type{T}, x::T) where {T<:INDICES} = x
Base.convert(::Type{Miller{T}}, mb::MillerBravais{T}) where {T<:RealSpace} =
    Miller{T}(2 * mb[1] + mb[2], 2 * mb[2] + mb[1], mb[4])
Base.convert(::Type{Miller{T}}, mb::MillerBravais{T}) where {T<:ReciprocalSpace} =
    Miller{T}(mb[1], mb[2], mb[4])
Base.convert(::Type{MillerBravais{T}}, m::Miller{T}) where {T<:RealSpace} =
    MillerBravais{T}(2 * m[1] - m[2], 2 * m[2] - m[1], -(m[1] + m[2]), 3 * m[3])
Base.convert(::Type{MillerBravais{T}}, m::Miller{T}) where {T<:ReciprocalSpace} =
    MillerBravais{T}(m[1], m[2], -(m[1] + m[2]), m[3])
