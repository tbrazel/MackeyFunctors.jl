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
    for (i, (K_idx, H_idx)) in enumerate(mf.context.covers)   # confirm which entry is H
        n   = subgroup_inclusion_index(mf.context, mf.context.covers[i])
        M   = mf.values[H_idx]
        res = mf.cover_restrictions[i]   # M(H) -> M(K)
        tr  = mf.cover_transfers[i]      # M(K) -> M(H)
        for g in gens(M)
            if tr(res(g)) != n * g
                return false
            end
        end
    end
    
    return true
end
