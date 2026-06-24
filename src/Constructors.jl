function constant_mackey_functor(context::MackeyContext, M::AbstractAlgebra.FPModule)
    R = base_ring(M)
    id_hom = identity_homomorphism(M)
    id_iso = identity_isomorphism(M)
    values = fill(M, length(context.subgroups))
    restrictions = fill(id_hom, length(context.covers))
    transfers = [R(subgroup_inclusion_index(context,cover_index))*id_hom for cover_index in context.covers]
    conjugations = fill(id_iso, length(context.generators), length(context.subgroups))
    MackeyFunctor(context, values, restrictions, transfers, conjugations)
end