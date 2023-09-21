using CrystallographyCore: AbstractCell
using Spglib:
    get_symmetry, get_ir_reciprocal_mesh, get_stabilized_reciprocal_mesh, eachpoint

export getsymmetry,
    get_reciprocal_mesh,
    get_stabilized_reciprocal_mesh,
    eachpoint,
    getstabilizer,
    getlittlegroup,
    getorbit,
    getstar

function getsymmetry(cell::AbstractCell, symprec=1e-5)
    rotations, translations = get_symmetry(cell, symprec)
    return map(rotations, translations) do rotation, translation
        SeitzOperator(rotation, translation)
    end
end

const get_reciprocal_mesh = get_ir_reciprocal_mesh

function getstabilizer(symops, point)
    return filter(symops) do symop
        isinvariant(symop, point)
    end
end
const getlittlegroup = getstabilizer

function getorbit(symops, point)
    return unique(
        map(symops) do symop
            symop(point)
        end,
    )
end
const getstar = getorbit

function isinvariant(op::SeitzOperator, 𝐫)
    𝐫′ = op(𝐫)
    if 𝐫′ == 𝐫
        return true
    elseif all(isinteger.(𝐫′ - 𝐫))
        return true
    else
        return false
    end
end
