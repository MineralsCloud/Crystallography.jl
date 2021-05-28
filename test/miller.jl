@testset "Test constructors" begin
    @test_throws TypeError Miller{Int,Int}
    @test_throws TypeError MillerBravais{Int,Int}
    @test_throws TypeError Miller{Int}(1, 2, 3)
    @test_throws TypeError MillerBravais{Int}(1, 2, 3, 4)
    @test_throws TypeError Miller{Int}([1, 2, 3])
    @test_throws TypeError MillerBravais{Int}([1, 2, 3, 4])
    @test_throws TypeError Miller{Int}((1, 2, 3))
    @test_throws TypeError MillerBravais{Int}((1, 2, 3, 4))
end

@testset "Test conversion between real `Miller` and `MillerBravais`" begin
    millerreal = [[1, 1, 1], [-1, 0, 1], [0, -1, 1], [-1, -1, 1], [1, 0, 1], [0, 1, 1]]
    millerbravaisreal = [
        [1, 1, -2, 3],
        [-2, 1, 1, 3],
        [1, -2, 1, 3],
        [-1, -1, 2, 3],
        [2, -1, -1, 3],
        [-1, 2, -1, 3],
    ]
    for (x, y) in zip(millerreal, millerbravaisreal)
        m = Miller{RealSpace}(x)
        mb = MillerBravais{RealSpace}(y)
        @test convert(MillerBravais{RealSpace}, m) == mb
        @test convert(Miller{RealSpace}, mb) == m
        @test convert(Miller{RealSpace}, m) == m
        @test convert(MillerBravais{RealSpace}, mb) == mb
        @test_throws MethodError convert(Miller{ReciprocalSpace}, m)
        @test_throws MethodError convert(MillerBravais{ReciprocalSpace}, mb)
        @test_throws MethodError convert(MillerBravais{ReciprocalSpace}, m)
        @test_throws MethodError convert(Miller{ReciprocalSpace}, mb)
        @test_throws MethodError convert(MillerBravais, m)
        @test_throws MethodError convert(Miller, mb)
    end
end

@testset "Test conversion between reciprocal `Miller` and `MillerBravais`" begin
    millerreciprocal =
        [[1, 0, 0], [0, 1, 0], [1, -1, 0], [-1, 0, 0], [0, -1, 0], [-1, 1, 0]]
    millerbravaisreciprocal = [
        [1, 0, -1, 0],
        [0, 1, -1, 0],
        [1, -1, 0, 0],
        [-1, 0, 1, 0],
        [0, -1, 1, 0],
        [-1, 1, 0, 0],
    ]
    for (x, y) in zip(millerreciprocal, millerbravaisreciprocal)
        m = Miller{ReciprocalSpace}(x)
        mb = MillerBravais{ReciprocalSpace}(y)
        @test convert(MillerBravais{ReciprocalSpace}, m) == mb
        @test convert(Miller{ReciprocalSpace}, mb) == m
        @test convert(Miller{ReciprocalSpace}, m) == m
        @test convert(MillerBravais{ReciprocalSpace}, mb) == mb
        @test_throws MethodError convert(Miller{RealSpace}, m)
        @test_throws MethodError convert(MillerBravais{RealSpace}, mb)
        @test_throws MethodError convert(MillerBravais{RealSpace}, m)
        @test_throws MethodError convert(Miller{RealSpace}, mb)
        @test_throws MethodError convert(MillerBravais, m)
        @test_throws MethodError convert(Miller, mb)
    end
end
