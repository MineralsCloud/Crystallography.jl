using LinearAlgebra: I, diagm

@testset "Test `SeitzOperator` constructors" begin
    𝐑 = diagm([1, 1, 1])
    𝐭 = [1, 2, 3]
    @testset "Matrix constructor" begin
        op = SeitzOperator(𝐑)
        @test getpointsymmetry(op) == 𝐑
        @test gettranslation(op) == zeros(3)
    end
    @testset "Vector constructor" begin
        op = SeitzOperator(𝐭)
        @test getpointsymmetry(op) == I
        @test gettranslation(op) == 𝐭
    end
    @testset "Matrix and vector constructor" begin
        op = SeitzOperator(𝐑, 𝐭)
        @test getpointsymmetry(op) == 𝐑
        @test gettranslation(op) == 𝐭
    end
    @testset "Construct by multiplication" begin
        𝐑₁ = [0 0 1; 1 0 0; 0 1 0]
        𝐭₁ = [1, 2, 3]
        𝐑₂ = [0 -1 0; 1 0 0; 0 0 1]
        𝐭₂ = [3, 2, 1]
        op₁ = SeitzOperator(𝐑₁, 𝐭₁)
        op₂ = SeitzOperator(𝐑₂, 𝐭₂)
        @test op₁ * op₂ isa SeitzOperator
        @test op₂ * op₁ isa SeitzOperator
    end
end

@testset "Test `ispointsymmetry`" begin
    @test ispointsymmetry(one(SeitzOperator{Int}))
    @test !ispointsymmetry(SeitzOperator([1, 2, 3]))
end

@testset "Test `istranslation`" begin
    𝐭 = [1, 2, 3]
    @test istranslation(SeitzOperator(𝐭))
    @test istranslation(one(SeitzOperator{Int}))
end

@testset "Composition of operators" begin
    𝐑₁ = [0 0 1; 1 0 0; 0 1 0]
    𝐭₁ = [1, 2, 3]
    𝐑₂ = [0 -1 0; 1 0 0; 0 0 1]
    𝐭₂ = [3, 2, 1]
    op₁ = SeitzOperator(𝐑₁, 𝐭₁)
    op₂ = SeitzOperator(𝐑₂, 𝐭₂)
    𝐫 = [4, 5, 6]
    result₁ = (op₂ ∘ op₁)(𝐫)
    result₂ = op₂(op₁(𝐫))
    @test result₁ == result₂
    @test result₁ != (op₁ ∘ op₂)(𝐫)
    @testset "Equivalence to matrix multiplication" begin
        op = op₂ * op₁
        result = op([𝐫; 1])
        @test result[1:3] == result₁
        op′ = op₁ * op₂
        result′ = op′([𝐫; 1])
        @test result != result′
        @test result′[1:3] == (op₁ ∘ op₂)(𝐫)
        @test_throws DimensionMismatch op([1, 2, 3, 4, 5])
    end
end

@testset "Test `inv`" begin
    𝐑 = diagm([1, 1, 1])
    𝐭 = [1, 2, 3]
    @test inv(SeitzOperator(𝐭)) == SeitzOperator(-𝐭)
    op = SeitzOperator(𝐑, 𝐭)
    @test inv(op) == SeitzOperator(𝐑, -𝐑 \ 𝐭)
    @test inv(inv(op)) ≈ op
    𝐫 = rand(3)
    𝐫′ = op(𝐫)
    @test inv(op)(𝐫′) ≈ 𝐫
    @test inv(one(SeitzOperator{Int})) == one(SeitzOperator{Int})
end

@testset "Test `one` and `oneunit`" begin
    T = SeitzOperator{Float64}
    @test one(T) == I
    @test isone(one(T))
    @test oneunit(T) == I
    @test isone(oneunit(T))
end

@testset "Test `similar`" begin
    @test typeof(similar(SeitzOperator{Float64})) == SeitzOperator{Float64}
    @test typeof(similar(SeitzOperator{Float64}, (4, 4))) == SeitzOperator{Float64}
end

# Example from https://uh.edu/~chembi/lect10sg.PDF
@testset "Test application of `SeitzOperator`s" begin
    𝐑 = [
        1 -1 0
        1 0 0
        0 0 1
    ]
    𝐭 = [1, 0, 0]
    op = SeitzOperator(𝐑, 𝐭)
    x, y, z = rand(3)
    @test op([x, y, z]) == [x - y + 1, x, z]
end

@testset "Test `conjugacy`" begin
    @testset "Example from Katayun Barmak's lecture notes" begin
        𝐱 = [1 / 2, 1 / 2, 0]
        W𝐭 = SeitzOperator(𝐱)
        W₆ = SeitzOperator([
            0 -1 0
            1 1 0
            0 0 1
        ])  # 6-fold rotation
        @test SeitzOperator(W₆, 𝐱) ==
            conjugacy(W𝐭, W₆) ==
            [
                0 -1 0 1
                1 1 0 -1/2
                0 0 1 0
                0 0 0 1
            ]
    end
    @testset "Along direction [0 0 1]" begin
        C₃ = SeitzOperator([
            0 -1 0
            1 -1 0
            0 0 1
        ])  # 3-fold rotation
        C₃² = SeitzOperator([
            -1 1 0
            -1 0 0
            0 0 1
        ])
        @test conjugacy(C₃, C₃²) == C₃²
    end
    @testset "Along direction [1 1 0]" begin
        C₂ = SeitzOperator([
            0 1 0
            1 0 0
            0 0 -1
        ])  # 2-fold rotation
        m = SeitzOperator([
            0 -1 0
            -1 0 0
            0 0 1
        ])  # Mirror
        @test conjugacy(C₂, m) == m
    end
end
