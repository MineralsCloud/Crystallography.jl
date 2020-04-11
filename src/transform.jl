struct CrystalCoord{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

struct RealFromReciprocal
    basis::SMatrix{3,3}
end
struct ReciprocalFromReal
    basis::SMatrix{3,3}
end
struct CartesianFromCrystal
    basis::SMatrix{3,3}
end
struct CrystalFromCartesian
    basis::SMatrix{3,3}
end
struct CrystalFromCrystal end
for T in
    (:RealFromReciprocal, :ReciprocalFromReal, :CartesianFromCrystal, :CrystalFromCartesian)
    eval(quote
        $T(m::AbstractMatrix) = $T(SMatrix{3,3}(m))
        $T(lattice::Lattice) = $T(convert(Matrix{eltype(lattice)}, lattice))
    end)
end

Base.inv(x::Union{CrystalFromCartesian,CartesianFromCrystal}) = typeof(x)(inv(x.basis))

(t::CrystalFromCartesian)(v::AbstractVector) =
    CrystalCoord(convert(Matrix{eltype(t.basis)}, t.basis) * v)
(t::CartesianFromCrystal)(v::CrystalCoord) =
    SVector(convert(Matrix{eltype(t.basis)}, t.basis)' * v)
(t::CrystalFromCrystal)(v::CrystalCoord) =
    CrystalFromCartesian(t.to)(CartesianFromCrystal(t.from)(v))

StaticArrays.similar_type(
    ::Type{<:CrystalCoord},  # Do not delete the `<:`!
    ::Type{T},
    size::Size{(3,)},
) where {T} = CrystalCoord{T}
