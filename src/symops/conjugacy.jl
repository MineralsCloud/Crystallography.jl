struct ConjugacyOperation{T}
    op::SeitzOperator{T}
end
# Faster than the other implementation
(conj::ConjugacyOperation)(op::SeitzOperator) = conj.op * op * inv(conj.op)

conjugacy_class(op::SeitzOperator, ops) = Set(ConjugacyOperation(opᵢ)(op) for opᵢ in ops)

function partition(group)
    return Set(
        map(group) do op
            conjugacy_class(op, filter(!=(op), group))
        end,
    )
end
