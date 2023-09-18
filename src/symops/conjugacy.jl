# Faster than the other implementation
conjugate(op₁::SeitzOperator, op₂::SeitzOperator) = op₁ * op₂ * inv(op₁)
