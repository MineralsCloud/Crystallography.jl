struct ConjugacyOperation{T}
    op::SeitzOperator{T}
end
# Faster than the other implementation
(conj::ConjugacyOperation)(op::SeitzOperator) = conj.op * op * inv(conj.op)
