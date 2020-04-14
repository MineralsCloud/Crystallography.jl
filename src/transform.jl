struct RealFromReciprocal
    m::SMatrix{3,3}
end
struct ReciprocalFromReal
    m::SMatrix{3,3}
end
struct CartesianFromCrystal
    m::SMatrix{3,3}
end
struct CrystalFromCartesian
    m::SMatrix{3,3}
end
for T in (:CartesianFromCrystal, :CrystalFromCartesian)
    eval(quote
        $T(lattice::Lattice) = $T(convert(Matrix{eltype(lattice)}, lattice))
    end)
end

Base.inv(x::CrystalFromCartesian) = CartesianFromCrystal(inv(x.m))
Base.inv(x::CartesianFromCrystal) = CrystalFromCartesian(inv(x.m))
Base.:∘(x::CrystalFromCartesian, y::CartesianFromCrystal) = x.m * y.m == I ? identity : error("undefined!")

(x::CartesianFromCrystal)(v::AbstractVector) = x.m * v
(x::CrystalFromCartesian)(v::AbstractVector) = inv(x.m) * v
