export Miller, MillerBravais, ReciprocalMiller, ReciprocalMillerBravais
export @m_str, @mb_str

abstract type Indices <: AbstractVector{Int} end
abstract type AbstractMiller <: Indices end
struct Miller <: AbstractMiller
    data::SVector{3,Int}
    Miller(v) = new(iszero(v) ? v : v .÷ gcd(v))
end
Miller(i, j, k) = Miller([i, j, k])
struct ReciprocalMiller <: AbstractMiller
    data::SVector{3,Int}
    ReciprocalMiller(v) = new(iszero(v) ? v : v .÷ gcd(v))
end
ReciprocalMiller(i, j, k) = ReciprocalMiller([i, j, k])

abstract type AbstractMillerBravais <: Indices end
struct MillerBravais <: AbstractMillerBravais
    data::SVector{4,Int}
    function MillerBravais(v)
        @assert v[3] == -v[1] - v[2] "the 3rd index of `MillerBravais` should equal to the negation of the first two!"
        return new(iszero(v) ? v : v .÷ gcd(v))
    end
end
MillerBravais(i, j, k, l) = MillerBravais([i, j, k, l])
struct ReciprocalMillerBravais <: AbstractMillerBravais
    data::SVector{4,Int}
    function ReciprocalMillerBravais(v)
        @assert v[3] == -v[1] - v[2] "the 3rd index of `MillerBravais` should equal to the negation of the first two!"
        return new(iszero(v) ? v : v .÷ gcd(v))
    end
end
ReciprocalMillerBravais(i, j, k, l) = ReciprocalMillerBravais([i, j, k, l])

# This is a helper function and should not be exported!
function _indices_str(r::Regex, s::AbstractString)
    m = match(r, strip(s))
    if m === nothing
        throw(ArgumentError("not a valid expression!"))
    else
        brackets = first(m.captures) * last(m.captures)
        x = (parse(Int, x) for x in m.captures[2:(end-1)])
        if brackets == "()"
            ReciprocalMiller(x...)
        elseif brackets == "{}"
            ReciprocalMillerBravais(x...)
        elseif brackets == "[]"
            Miller(x...)
        elseif brackets == "<>"
            MillerBravais(x...)
        else
            @assert false "this should never happen!"
        end
    end
end

macro m_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    return _indices_str(r, s)
end

macro mb_str(s)
    r =
        r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    return _indices_str(r, s)
end

Base.size(::AbstractMiller) = (3,)
Base.size(::AbstractMillerBravais) = (4,)

Base.IndexStyle(::Type{<:Indices}) = IndexLinear()

Base.getindex(x::Indices, i) = getindex(x.data, i)

Base.convert(::Type{T}, x::T) where {T<:Indices} = x
Base.convert(::Type{Miller}, mb::MillerBravais) =
    Miller(2 * mb[1] + mb[2], 2 * mb[2] + mb[1], mb[4])
Base.convert(::Type{ReciprocalMiller}, mb::ReciprocalMillerBravais) =
    ReciprocalMiller(mb[1], mb[2], mb[4])
Base.convert(::Type{MillerBravais}, m::Miller) =
    MillerBravais(2 * m[1] - m[2], 2 * m[2] - m[1], -(m[1] + m[2]), 3 * m[3])
Base.convert(::Type{ReciprocalMillerBravais}, m::ReciprocalMiller) =
    ReciprocalMillerBravais(m[1], m[2], -(m[1] + m[2]), m[3])
