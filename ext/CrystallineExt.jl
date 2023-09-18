module CrystallineExt

using Crystallography: getpointsymmetry, gettranslation
using Crystalline: seitz

import Crystalline: SymOperation
import Crystallography: SeitzOperator

SeitzOperator(op::SymOperation{3}) = SeitzOperator(op.rotation, op.translation)

SymOperation(op::SeitzOperator) = SymOperation{3}(getpointsymmetry(op), gettranslation(op))

function Base.show(io::IO, op::SeitzOperator)
    if get(io, :compact, false)
        print(io, seitz(SymOperation(op)))
    else
        print(io, "SeitzOperator(", parent(op), ')')
    end
end
function Base.show(io::IO, ::MIME"text/plain", op::SeitzOperator)
    println(io, string(typeof(op)), " (", seitz(SymOperation(op)), ')')
    for row in eachrow(parent(op))
        println(io, ' ', join(row, "  "))
    end
end

end
