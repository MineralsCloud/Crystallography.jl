using CrystallographyCore: AbstractCell
using Spglib: get_symmetry

export list_symmetry

function list_symmetry(cell::AbstractCell, symprec=1e-5)
    rotations, translations = get_symmetry(cell, symprec)
    return map(rotations, translations) do rotation, translation
        SeitzOperator(rotation, translation)
    end
end
