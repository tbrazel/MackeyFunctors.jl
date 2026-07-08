using GAP
using AbstractAlgebra
using MackeyFunctors



"""
The following method is a boolean to check if a Mackey functor is cohomological.  
A Mackey functor M is cohomological if for every cover relation K ≤ H,
the composition of transfer followed by restriction equals multiplication
by the index [H : K] on M(K).  
"""
function is_cohomological(mf::MackeyFunctor)
    for (i, (K_idx, _)) in enumerate(mf.context.covers)
        inclusion_index = subgroup_inclusion_index(mf.context, mf.context.covers[i])

        lhs = mf.cover_restrictions[i] * mf.cover_transfers[i]
        rhs = inclusion_index * identity_map(mf.values[K_idx])

        if lhs != rhs
            return false
        end
    end

    return true
end
