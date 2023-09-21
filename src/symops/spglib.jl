using CrystallographyCore: AbstractCell
using Spglib: get_symmetry

export findsymmetry, findstabilizer, findlittlegroup, findorbit, findstar

function findsymmetry(cell::AbstractCell, symprec=1e-5)
    rotations, translations = get_symmetry(cell, symprec)
    return map(rotations, translations) do rotation, translation
        SeitzOperator(rotation, translation)
    end
end

function findstabilizer(symops, point, rtol=eps())
    return filter(symops) do symop
        isapprox(symop(point), point; rtol=rtol)
    end
end
const findlittlegroup = findstabilizer

function findorbit(symops, point)
    return unique(
        map(symops) do symop
            symop(point)
        end,
    )
end
const findstar = findorbit
