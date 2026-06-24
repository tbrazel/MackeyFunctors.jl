function constant_mackey_functor(context::MackeyContext, M::AbstractAlgebra.FPModule)
    R = base_ring(M)
    id_hom = identity_homomorphism(M)
    id_iso = identity_isomorphism(M)
    values = fill(M, length(context.subgroups))
    restrictions = fill(id_hom, length(context.covers))
    transfers = map(context.covers) do (i, j)
        H = context.subgroups[i]
        K = context.subgroups[j]
        R(GAP.Globals.Index(K, H)) * id_hom
    end
    conjugations = fill(id_iso, length(context.generators), length(context.subgroups))
    MackeyFunctor(context, values, restrictions, transfers, conjugations)
end