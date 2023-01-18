export FractionalFromCartesian,
    CartesianFromFractional,
    FractionalToCartesian,
    CartesianToFractional,
    StandardizedFromPrimitive,
    PrimitiveFromStandardized,
    PrimitiveToStandardized,
    StandardizedToPrimitive

abstract type ChangeOfBasis{T} <: AbstractMatrix{T} end
struct CartesianFromFractional{T} <: ChangeOfBasis{T}
    tf::SMatrix{3,3,T,9}
end
struct FractionalFromCartesian{T} <: ChangeOfBasis{T}
    tf::SMatrix{3,3,T,9}
end
# This requires the a-vector is parallel to the Cartesian x-axis.
# See https://en.wikipedia.org/wiki/Fractional_coordinates
"""
    CartesianFromFractional(lattice::Union{Lattice,ReciprocalLattice})
    CartesianFromFractional(a, b, c, α, β, γ)

Get the transformation from fractional coordinates to Cartesian coordinates.
"""
CartesianFromFractional(lattice::Lattice) = CartesianFromFractional(lattice.data)
CartesianFromFractional(lattice::ReciprocalLattice) =
    CartesianFromFractional(transpose(lattice.data))
function CartesianFromFractional(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    b_sinγ, b_cosγ = b .* (sind(γ), cosd(γ))
    return CartesianFromFractional(
        [
            a b_cosγ c*cosd(β)
            0 b_sinγ c*_auxiliary(α, β, γ)
            0 0 Ω/(a * b_sinγ)
        ]
    )
end
"""
    FractionalFromCartesian(lattice::Union{Lattice,ReciprocalLattice})
    FractionalFromCartesian(a, b, c, α, β, γ)

Get the transformation from Cartesian coordinates to fractional coordinates.
"""
FractionalFromCartesian(lattice::Lattice) = FractionalFromCartesian(inv(lattice.data))
FractionalFromCartesian(lattice::ReciprocalLattice) =
    FractionalFromCartesian(transpose(inv(lattice.data)))
function FractionalFromCartesian(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    b_sinγ = b * sind(γ)
    return FractionalFromCartesian(
        [
            1/a -cotd(γ)/a -b * c * _auxiliary(β, α, γ)/Ω
            0 1/b_sinγ -a * c * _auxiliary(α, β, γ)/Ω
            0 0 a * b_sinγ/Ω
        ],
    )
end
const FractionalToCartesian = CartesianFromFractional
const CartesianToFractional = FractionalFromCartesian

# This is a helper function and should not be exported!
_auxiliary(α, β, γ) = (cosd(α) - cosd(β) * cosd(γ)) / sind(γ)

(x::Union{CartesianFromFractional,FractionalFromCartesian})(v) = x.tf * collect(v)
(x::Union{CartesianFromFractional,FractionalFromCartesian})(p::ReciprocalPoint) =
    ReciprocalPoint(x.tf * collect(p.coord), p.weight)

Base.inv(x::FractionalFromCartesian) = CartesianFromFractional(inv(x.tf))
Base.inv(x::CartesianFromFractional) = FractionalFromCartesian(inv(x.tf))
Base.:∘(x::CartesianFromFractional, y::FractionalFromCartesian) = ∘(y, x)
Base.:∘(x::FractionalFromCartesian, y::CartesianFromFractional) =
    x.tf * y.tf ≈ I ? identity : error("undefined!")

# Idea from https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell
struct StandardizedFromPrimitive{T} <: ChangeOfBasis{T}
    tf::SMatrix{3,3,T,9}
end
struct PrimitiveFromStandardized{T} <: ChangeOfBasis{T}
    tf::SMatrix{3,3,T,9}
end
"""
    PrimitiveFromStandardized(tf::AbstractMatrix)

Construct the transformation from a standardized cell to a primitive cell.
"""
PrimitiveFromStandardized(tf::AbstractMatrix) = PrimitiveFromStandardized{eltype(tf)}(tf)
"""
    StandardizedFromPrimitive(tf::AbstractMatrix)

Construct the transformation from a primitive cell to a standardized cell.
"""
StandardizedFromPrimitive(tf::AbstractMatrix) = StandardizedFromPrimitive{eltype(tf)}(tf)
const PrimitiveToStandardized = StandardizedFromPrimitive
const StandardizedToPrimitive = PrimitiveFromStandardized

(x::Union{StandardizedFromPrimitive,PrimitiveFromStandardized})(v) = inv(x.tf) * collect(v)
(x::Union{StandardizedFromPrimitive,PrimitiveFromStandardized})(lattice::Lattice) =
    Lattice(lattice.data * x.tf)

Base.inv(x::StandardizedFromPrimitive) = PrimitiveFromStandardized(inv(x.tf))
Base.inv(x::PrimitiveFromStandardized) = StandardizedFromPrimitive(inv(x.tf))
Base.:∘(x::PrimitiveFromStandardized, y::StandardizedFromPrimitive) = ∘(y, x)
Base.:∘(x::StandardizedFromPrimitive, y::PrimitiveFromStandardized) =
    x.tf * y.tf ≈ I ? identity : error("undefined!")

Base.size(::ChangeOfBasis) = (3, 3)

Base.getindex(x::ChangeOfBasis, i) = getindex(x.tf, i)

Base.IndexStyle(::Type{<:ChangeOfBasis}) = IndexLinear()

function Base.show(io::IO, x::ChangeOfBasis)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        for row in eachrow(x.tf)
            println(io, ' ', join(row, "  "))
        end
    end
end
