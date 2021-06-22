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
        m = Miller(x)
        mb = MillerBravais(y)
        @test convert(MillerBravais, m) == mb
        @test convert(Miller, mb) == m
        @test convert(Miller, m) == m
        @test convert(MillerBravais, mb) == mb
        @test_throws MethodError convert(ReciprocalMiller, m)
        @test_throws MethodError convert(ReciprocalMillerBravais, mb)
        @test_throws MethodError convert(ReciprocalMillerBravais, m)
        @test_throws MethodError convert(ReciprocalMiller, mb)
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
        m = ReciprocalMiller(x)
        mb = ReciprocalMillerBravais(y)
        @test convert(ReciprocalMillerBravais, m) == mb
        @test convert(ReciprocalMiller, mb) == m
        @test convert(ReciprocalMiller, m) == m
        @test convert(ReciprocalMillerBravais, mb) == mb
        @test_throws MethodError convert(Miller, m)
        @test_throws MethodError convert(MillerBravais, mb)
        @test_throws MethodError convert(MillerBravais, m)
        @test_throws MethodError convert(Miller, mb)
    end
end
