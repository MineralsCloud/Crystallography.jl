using PeriodicTable: elements
using Unitful: uconvert

# Data from https://materialsproject.org/materials/mp-2657/
@testset "Test `crystaldensity`" begin
    lattice =
        Lattice(4.653272u"angstrom", 4.653272u"angstrom", 2.969203u"angstrom", 90, 90, 90)
    positions = [
        [0, 0, 0],
        [0.5, 0.5, 0.5],
        [0.19542, 0.80458, 0.5],
        [0.80458, 0.19542, 0.5],
        [0.30458, 0.30458, 0],
        [0.69542, 0.69542, 0],
    ]
    @testset "Test `Element`s" begin
        atoms = elements[[22, 22, 8, 8, 8, 8]]
        cell = Cell(lattice, positions, atoms)
        @test cellvolume(cell) ≈ 64.29197531534862u"angstrom^3"
        @test uconvert(u"g/cm^3", crystaldensity(lattice, atoms)) ≈
              4.125526333805304u"g/cm^3"
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.125526333805304u"g/cm^3"
    end
    @testset "Test `AbstractString`s" begin
        atoms = elements[["Titanium", "Titanium", "Oxygen", "Oxygen", "Oxygen", "Oxygen"]]
        cell = Cell(lattice, positions, atoms)
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.125526333805304u"g/cm^3"
    end
    @testset "Test `Integer`s" begin
        atoms = [22, 22, 8, 8, 8, 8]
        cell = Cell(lattice, positions, atoms)
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.125526333805304u"g/cm^3"
    end
    @testset "Test `Symbol`s" begin
        atoms = [:Ti, :Ti, :O, :O, :O, :O]
        cell = Cell(lattice, positions, atoms)
        @test uconvert(u"g/cm^3", crystaldensity(cell)) ≈ 4.125526333805304u"g/cm^3"
    end
end
