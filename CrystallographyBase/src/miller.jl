using Combinatorics: permutations

export Miller, MillerBravais, ReciprocalMiller, ReciprocalMillerBravais
export family, @m_str

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
        if brackets in ("()", "{}")
            if length(x) == 3
                ReciprocalMiller(x...)
            else  # length(x) == 4
                ReciprocalMillerBravais(x...)
            end
        elseif brackets ∈ ("[]", "<>")
            if length(x) == 3
                Miller(x...)
            else  # length(x) == 4
                MillerBravais(x...)
            end
        else
            @assert false "this should never happen!"
        end
    end
end

macro m_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)?[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    return _indices_str(r, s)
end

function family(mb::T) where {T<:AbstractMillerBravais}
    perm = collect(permutations(mb[1:3]))  # Permute the first 3 indices for equivalent basis vectors
    negate = -perm  # Use negative indices
    pool = unique(append!(perm, negate))
    allowed = filter(v -> v[3] == -v[1] - v[2], pool)
    return map(allowed) do x
        T(x..., mb[4])  # Add the 4th index back
    end
end
function family(m::T) where {T<:AbstractMiller}
    mb = convert(T <: Miller ? MillerBravais : ReciprocalMillerBravais, m)  # Real or reciprocal space
    vec = family(mb)
    return map(x -> convert(T, x), vec)
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

function Base.show(io::IO, x::Union{Miller,MillerBravais})
    print(io, '<')
    print(io, join(x.data, " "))
    print(io, '>')
end
function Base.show(io::IO, x::Union{ReciprocalMiller,ReciprocalMillerBravais})
    print(io, '{')
    print(io, join(x.data, " "))
    print(io, '}')
end
