# Examples

```@contents
Pages = ["examples.md"]
Depth = 2
```

## Creating a `Cell`

To create a `Cell`, we first need to create a `Lattice`.
There are multiple ways of doing it. For example, if we know the six lattice constants,
we can do

```@repl 1
using CrystallographyBase, Unitful, UnitfulAtomic
lattice₁ = Lattice(4u"nm", 180u"bohr", 3u"angstrom", 90, 90, 90)
```

Or, equivalently,

```@repl 1
lattice₁ == Lattice([
    4u"nm" 0u"m" 0.0u"cm"
    0u"cm" 180.0u"bohr" 0u"m"
    0u"bohr" 0u"nm" (3//1)*u"angstrom"
])
```

Then we can add atoms and their positions (in crystal coordinates):

```@repl 1
lattice₂ = [
    -3.0179389205999998 -3.0179389205999998 0.0000000000000000
    -5.2272235447000002 5.2272235447000002 0.0000000000000000
    0.0000000000000000 0.0000000000000000 -9.7736219469000005
]
positions = [[2 / 3, 1 / 3, 1 / 4], [1 / 3, 2 / 3, 3 / 4]]
atoms = [1, 1]
cell = Cell(lattice₂, positions, atoms)
```

## Reciprocal space

To get the reciprocal lattice, we run `reciprocal`:

```@repl 1
reciprocal(lattice₁)
```

!!! note
    Never use `ReciprocalLattice` directly unless you know what you are doing!

## Supercell generation

We can specify the replication factors in each direction in the following ways:

```@repl 1
supercell(lattice₁, [2, 3, 4])
supercell(lattice₁, 3)
supercell(cell, [2, 3, 4])
supercell(cell, 3)
```

If only one integer is provided, it will be used in all three spatial directions.
