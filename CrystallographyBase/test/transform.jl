using CoordinateTransformations: IdentityTransformation

@testset "Fractional coordinates to Cartesian coordinates" begin
    lattice = Lattice([
        1/2 0 0
        0 1/2 0
        0 0 1
    ])
    @test CartesianFromFractional(lattice)([2, 3, 1]) == [1, 3 / 2, 1]
    # Compared with "SeeK-path" example oC1 (SiTi)
    @testset "Base centered orthorhombic" begin
        a, b, c = 3.67939899, 4.15357986, 3.8236792350
        lattice = Lattice([a, -b, 0] / 2, [a, b, 0] / 2, [0, 0, c])
        f2c = CartesianFromFractional(lattice)
        @test f2c([1 / 2, 1 / 2, 1 / 2]) == [1.8396994950, 0, 1.9118396175]  # Si
        @test f2c([0, 0, 0]) == [0, 0, 0]  # Ti
        c2f = inv(f2c)
        @test c2f([1.8396994950, 0, 1.9118396175]) ≈ [1 / 2, 1 / 2, 1 / 2]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
end

# Example from http://ww1.iucr.org/iucr-top/comm/cteach/pamphlets/22/node34.html#SECTION00073300000000000000
@testset "Refer ZrSiO₄ primitive cell with origin at 2/m" begin
    𝐚ₚ, 𝐛ₚ, 𝐜ₚ = [1, 0, 0], [0, 1, 0], [1, 1, 1] // 2
    P = PrimitiveFromStandardized(hcat(𝐚ₚ, 𝐛ₚ, 𝐜ₚ))
    @test P([0, 1 // 4, 7 // 8]) == [-7, -5, 14] // 8  # The new coordinates of the first Zr atom
    @test P([0, 1 // 4, 3 // 8]) == [-3, -1, 6] // 8  # The new coordinates of the first Si atom
    @test P([0, 0.45, 0.215]) ≈ [-0.215, 0.235, 0.43]  # The new coordinates of the first O atom
    P⁻¹ = inv(P)
    @test P⁻¹([-7, -5, 14] // 8) == [0, 1 // 4, 7 // 8]
    @test P⁻¹([-3, -1, 6] // 8) == [0, 1 // 4, 3 // 8]
    @test P⁻¹([-0.215, 0.235, 0.43]) ≈ [0, 0.45, 0.215]
    @test P ∘ P⁻¹ == P⁻¹ ∘ P == IdentityTransformation()
end

@testset "Special points in the Brillouin zone" begin  # Reciprocal space transformation
    # See http://lampx.tugraz.at/~hadley/ss1/bzones/orthorhombic.php
    @testset "Simple orthorhombic Brillouin zone" begin
        a, b, c = 2, 3, 5
        reci_lattice = inv(Lattice(a, b, c, 90, 90, 90))
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([1 / 2, 0, 0]) == [1 / a, 0, 0] / 2
        @test f2c([0, 1 / 2, 0]) == [0, 1 / b, 0] / 2
        @test f2c([0, 0, 1 / 2]) == [0, 0, 1 / c] / 2
        @test f2c([0, 1 / 2, 1 / 2]) == [0, 1 / b, 1 / c] / 2
        @test f2c([1 / 2, 1 / 2, 0]) == [1 / a, 1 / b, 0] / 2
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([1 / a, 0, 0] / 2) == [1 / 2, 0, 0]
        @test c2f([0, 1 / b, 0] / 2) == [0, 1 / 2, 0]
        @test c2f([0, 0, 1 / c] / 2) == [0, 0, 1 / 2]
        @test c2f([0, 1 / b, 1 / c] / 2) == [0, 1 / 2, 1 / 2]
        @test c2f([1 / a, 1 / b, 0] / 2) == [1 / 2, 1 / 2, 0]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
    # Compared with "SeeK-path" example oC1 (SiTi)
    # The results on http://lampx.tugraz.at/~hadley/ss1/bzones/orthorhombic_bc.php are wrong by a factor of 2.
    @testset "Base centered orthorhombic Brillouin zone" begin
        a, b, c = 3.67939899, 4.15357986, 3.8236792350
        reci_lattice = inv(Lattice([a, -b, 0] / 2, [a, b, 0] / 2, [0, 0, c]))
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([-1 / 2, 1 / 2, 0]) ≈ [0, 1 / b, 0] ≈ [0, 1.5127156619, 0] / 2pi  # Y'
        @test f2c([0, 0, 1 / 2]) ≈ [0, 0, 1 / 2 / c] ≈ [0, 0, 0.8216151148] / 2pi  # Z
        @test f2c([-1 / 2, 1 / 2, 1 / 2]) ≈
              [0, 1 / b, 1 / 2 / c] ≈
              [0, 1.5127156619, 0.8216151148] / 2pi  # T'
        @test f2c([0, 1 / 2, 0]) ≈
              [1 / 2 / a, 1 / 2 / b, 0] ≈
              [0.8538331021, 0.7563578310, 0] / 2pi  # S
        @test f2c([0, 1 / 2, 1 / 2]) ≈
              [1 / a, 1 / b, 1 / c] / 2 ≈
              [0.8538331021, 0.7563578310, 0.8216151148] / 2pi # R
        @test f2c([-0.4461772527, 0.5538227473, 1 / 2]) ≈
              [0.1838225731, 1.5127156619, 0.8216151148] / 2pi  # E0
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([0, 1 / b, 0]) ≈ [-1 / 2, 1 / 2, 0]
        @test c2f([0, 0, 1 / 2 / c]) ≈ [0, 0, 1 / 2]
        @test c2f([0, 1 / b, 1 / 2 / c]) ≈ [-1 / 2, 1 / 2, 1 / 2]
        @test c2f([1 / 2 / a, 1 / 2 / b, 0]) ≈ [0, 1 / 2, 0]
        @test c2f([1 / a, 1 / b, 1 / c] / 2) ≈ [0, 1 / 2, 1 / 2]
        @test c2f([0.1838225731, 1.5127156619, 0.8216151148] / 2pi) ≈
              [-0.4461772527, 0.5538227473, 1 / 2]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
    # See http://lampx.tugraz.at/~hadley/ss1/bzones/tetbc.php
    @testset "Body centered tetragonal Brillouin zone" begin
        a, c = 4, 6
        reci_lattice = inv(Lattice([a, a, -c] / 2, [a, -a, c] / 2, [-a, a, c] / 2))
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([0, 0, 0]) == [0, 0, 0]
        @test f2c([1 / 2, 0, 0]) == [1 / a, 1 / a, 0] / 2
        @test f2c([1 / 2, 1 / 2, -1 / 2]) == [1 / a, 0, 0]
        @test f2c([0, 1 / 2, 0]) == [1 / a, 0, 1 / c] / 2
        @test f2c([1, 1, 1] / 4) == [1 / a, 1 / a, 1 / c] / 2
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([0, 0, 0]) == [0, 0, 0]
        @test c2f([1 / a, 1 / a, 0] / 2) == [1 / 2, 0, 0]
        @test c2f([1 / a, 0, 0]) == [1 / 2, 1 / 2, -1 / 2]
        @test c2f([1 / a, 0, 1 / c] / 2) ≈ [0, 1 / 2, 0]
        @test c2f([1 / a, 1 / a, 1 / c] / 2) == [1, 1, 1] / 4
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
    # See http://lampx.tugraz.at/~hadley/ss1/bzones/hexagonal.php
    @testset "Simple hexagonal Brillouin zone" begin
        a, c = 2, 3.2
        reci_lattice = inv(Lattice([a, 0, 0], [a / 2, sqrt(3) / 2 * a, 0], [0, 0, c]))
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([0, 0, 0]) == [0, 0, 0]  # Γ
        @test f2c([1 / 2, 0, 0]) ≈ [1 / a, -1 / sqrt(3) / a, 0] / 2  # M
        @test f2c([0, 0, 1 / 2]) ≈ [0, 0, 1 / c] / 2  # A
        @test f2c([2 / 3, 1 / 3, 0]) ≈ [2 / 3 / a, 0, 0]  # K
        @test f2c([2 / 3, 1 / 3, 1 / 2]) ≈ [2 / 3 / a, 0, 1 / 2 / c]  # H
        @test f2c([1 / 2, 0, 1 / 2]) ≈ [1 / a, -1 / sqrt(3) / a, 1 / c] / 2  # L
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([0, 0, 0]) == [0, 0, 0]
        @test c2f([1 / a, -1 / sqrt(3) / a, 0] / 2) == [1 / 2, 0, 0]
        @test c2f([0, 0, 1 / c] / 2) == [0, 0, 1 / 2]
        @test c2f([2 / 3 / a, 0, 0]) == [2 / 3, 1 / 3, 0]
        @test c2f([2 / 3 / a, 0, 1 / 2 / c]) == [2 / 3, 1 / 3, 1 / 2]
        @test c2f([1 / a, -1 / sqrt(3) / a, 1 / c] / 2) == [1 / 2, 0, 1 / 2]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
    # See http://lampx.tugraz.at/~hadley/ss1/bzones/sc.php
    @testset "Simple cubic Brillouin zone" begin
        a = 4
        reci_lattice = inv(Lattice(a, a, a, 90, 90, 90))
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([0, 0, 0]) == [0, 0, 0]
        @test f2c([1 / 2, 1 / 2, 0]) == [1 / a, 1 / a, 0] / 2
        @test f2c([1 / 2, 1 / 2, 1 / 2]) == [1 / a, 1 / a, 1 / a] / 2
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([0, 0, 0]) == [0, 0, 0]
        @test c2f([1 / a, 1 / a, 0] / 2) == [1 / 2, 1 / 2, 0]
        @test c2f([1 / a, 1 / a, 1 / a] / 2) == [1 / 2, 1 / 2, 1 / 2]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
    # See http://lampx.tugraz.at/~hadley/ss1/bzones/fcc.php
    @testset "Face centered cubic Brillouin zone" begin
        a = 4
        reci_lattice = inv(Lattice([
            1 1 0
            0 1 1
            1 0 1
        ]) * a / 2)
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([0, 0, 0]) == [0, 0, 0]
        @test f2c([0, 1 / 2, 1 / 2]) == [0, 1 / a, 0]
        @test f2c([1 / 2, 1 / 2, 1 / 2]) == [1 / a, 1 / a, 1 / a] / 2
        @test f2c([1 / 4, 3 / 4, 1 / 2]) == [1 / 2 / a, 1 / a, 0]
        @test f2c([1 / 4, 5 / 8, 5 / 8]) == [1 / 2 / a, 2 / a, 1 / 2 / a] / 2
        @test f2c([3 / 8, 3 / 4, 3 / 8]) == [3 / 2 / a, 3 / 2 / a, 0] / 2
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([0, 0, 0]) == [0, 0, 0]
        @test f2c([0, 1 / 2, 1 / 2]) == [0, 1 / a, 0]
        @test f2c([1 / 2, 1 / 2, 1 / 2]) == [1 / a, 1 / a, 1 / a] / 2
        @test c2f([1 / 2 / a, 1 / a, 0]) == [1 / 4, 3 / 4, 1 / 2]
        @test c2f([1 / 2 / a, 2 / a, 1 / 2 / a] / 2) == [1 / 4, 5 / 8, 5 / 8]
        @test c2f([3 / 2 / a, 3 / 2 / a, 0] / 2) == [3 / 8, 3 / 4, 3 / 8]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
    # See http://lampx.tugraz.at/~hadley/ss1/bzones/bcc.php
    @testset "Body centered cubic Brillouin zone" begin
        a = 4
        reci_lattice = inv(Lattice([
            1 -1 1
            1 1 -1
            -1 1 1
        ]) * a / 2)
        f2c = CartesianFromFractional(reci_lattice)
        @test f2c([0, 0, 0]) == [0, 0, 0]
        @test f2c([-1 / 2, 1 / 2, 1 / 2]) == [0, 0, 1 / a]
        @test f2c([1 / 4, 1 / 4, 1 / 4]) == [1 / a, 1 / a, 1 / a] / 2
        @test f2c([0, 1 / 2, 0]) == [0, 1 / a, 1 / a] / 2
        c2f = FractionalFromCartesian(reci_lattice)
        @test c2f([0, 0, 0]) == [0, 0, 0]
        @test c2f([0, 0, 1 / a]) == [-1 / 2, 1 / 2, 1 / 2]
        @test c2f([1 / a, 1 / a, 1 / a] / 2) == [1 / 4, 1 / 4, 1 / 4]
        @test c2f([0, 1 / a, 1 / a] / 2) == [0, 1 / 2, 0]
        @test c2f ∘ f2c == f2c ∘ c2f == IdentityTransformation()
    end
end
