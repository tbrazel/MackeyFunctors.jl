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
    subgroups = mf.context.subgroups

    for (K_idx, H_idx) in mf.context.covers
        K = subgroups[K_idx]
        H = subgroups[H_idx]

        lhs = mf.cover_restrictions[(H_idx, K_idx)] * mf.cover_transfers[(K_idx, H_idx)]  # ✅ use mf
        rhs = index(H, K) * identity_map(mf.modules[K_idx])

        if lhs != rhs
            return false
        end
    end

    return true
end
