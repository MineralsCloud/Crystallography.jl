using CoordinateTransformations: IdentityTransformation, Transformation

export FractionalFromCartesian,
    CartesianFromFractional,
    FractionalToCartesian,
    CartesianToFractional,
    StandardizedFromPrimitive,
    PrimitiveFromStandardized,
    PrimitiveToStandardized,
    StandardizedToPrimitive

struct CartesianFromFractional{T} <: Transformation
    tf::SMatrix{3,3,T,9}
end
struct FractionalFromCartesian{T} <: Transformation
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
    x.tf * y.tf ≈ I ? IdentityTransformation() : error("undefined!")

# Idea from https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell
struct StandardizedFromPrimitive{T} <: Transformation
    tf::SMatrix{3,3,T,9}
end
struct PrimitiveFromStandardized{T} <: Transformation
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
(x::PrimitiveFromStandardized)(lattice::Lattice) = Lattice(lattice.data * x.tf)
(x::StandardizedFromPrimitive)(lattice::Lattice) = Lattice(lattice.data * inv(x.tf))

# PrimitiveFromStandardized(::ACentering) = PrimitiveFromStandardized([
#     1 0 0
#     0 1//2 -1//2
#     0 1//2 1//2
# ])
# PrimitiveFromStandardized(::CCentering) = PrimitiveFromStandardized([
#     1//2 1//2 0
#     -1//2 1//2 0
#     0 0 1
# ])
# PrimitiveFromStandardized(::RhombohedralCentering) = PrimitiveFromStandardized([
#     2//3 -1//3 -1//3
#     1//3 1//3 -2//3
#     1//3 1//3 1//3
# ])
# PrimitiveFromStandardized(::BodyCentering) = PrimitiveFromStandardized([
#     -1//2 1//2 1//2
#     1//2 -1//2 1//2
#     1//2 1//2 -1//2
# ])
# PrimitiveFromStandardized(::FaceCentering) = PrimitiveFromStandardized([
#     0 1//2 1//2
#     1//2 0 1//2
#     1//2 1//2 0
# ])
# StandardizedFromPrimitive(::ACentering) = StandardizedFromPrimitive([
#     1 0 0
#     0 1 1
#     0 -1 1
# ])
# StandardizedFromPrimitive(::CCentering) = StandardizedFromPrimitive([
#     1 -1 0
#     1 1 0
#     0 0 1
# ])
# StandardizedFromPrimitive(::RhombohedralCentering) = StandardizedFromPrimitive([
#     1 0 1
#     -1 1 1
#     0 -1 1
# ])
# StandardizedFromPrimitive(::BodyCentering) = StandardizedFromPrimitive([
#     0 1 1
#     1 0 1
#     1 1 0
# ])
# StandardizedFromPrimitive(::FaceCentering) = StandardizedFromPrimitive([
#     -1 1 1
#     1 -1 1
#     1 1 -1
# ])

Base.inv(x::StandardizedFromPrimitive) = PrimitiveFromStandardized(inv(x.tf))
Base.inv(x::PrimitiveFromStandardized) = StandardizedFromPrimitive(inv(x.tf))
Base.:∘(x::PrimitiveFromStandardized, y::StandardizedFromPrimitive) = ∘(y, x)
Base.:∘(x::StandardizedFromPrimitive, y::PrimitiveFromStandardized) =
    x.tf * y.tf ≈ I ? IdentityTransformation() : error("undefined!")

const ChangeOfBasis{T} = Union{
    CartesianFromFractional{T},
    FractionalFromCartesian{T},
    StandardizedFromPrimitive{T},
    PrimitiveFromStandardized{T},
}
Base.iterate(x::ChangeOfBasis) = iterate(x.tf)
Base.iterate(x::ChangeOfBasis, state) = iterate(x.tf, state)

Base.eltype(::Type{<:ChangeOfBasis{T}}) where {T} = T

Base.length(::ChangeOfBasis) = 9

Base.size(::ChangeOfBasis) = (3, 3)
Base.size(::ChangeOfBasis, dim::Integer) = dim <= 2 ? 3 : 1

Base.IteratorSize(::Type{<:ChangeOfBasis}) = Base.HasShape{2}()

# Enables `firstindex(x, dim)` and `x[1, 2:end]` or `x[begin:2, 2]`
Base.axes(x::ChangeOfBasis, dim::Integer) = axes(x.tf, dim)

Base.getindex(x::ChangeOfBasis, i) = getindex(x.tf, i)
Base.getindex(x::ChangeOfBasis, I::Vararg) = getindex(x.tf, I...)

Base.firstindex(::ChangeOfBasis) = 1

Base.lastindex(::ChangeOfBasis) = 9

function Base.show(io::IO, x::ChangeOfBasis)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        for row in eachrow(x.tf)
            println(io, ' ', join(row, " "))
        end
    end
end
