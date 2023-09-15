module CrystallineExt

using Crystallography: getpointsymmetry, gettranslation

import Crystalline: SymOperation
import Crystallography: SeitzOperator

SeitzOperator(op::SymOperation{3}) = SeitzOperator(op.rotation, op.translation)

SymOperation(op::SeitzOperator) = SymOperation{3}(getpointsymmetry(op), gettranslation(op))

end
