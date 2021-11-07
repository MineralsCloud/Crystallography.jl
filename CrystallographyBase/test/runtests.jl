using CrystallographyBase
using Test

@testset "CrystallographyBase.jl" begin
    include("lattice.jl")
    include("miller.jl")
    include("metric.jl")
    include("transform.jl")
    include("reciprocal.jl")
    # include("wignerseitz.jl")
end
