"""
    constant_mackey_functor(ctx,M)
    constant_mackey_functor(ctx,R)

Given a `MackeyContext` and an ``R``-module ``M``, this method outputs the fixed-point Mackey functor for ``M`` with trivial ``G``-action. This is also called the *constant Mackey functor* valued at ``M``.

This method can also be fed a context and a ring ``R``, and it will output the constant Mackey functor valued at ``R`` considered as a free rank one ``R``-module.
"""
function constant_mackey_functor end

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

function constant_mackey_functor(context::MackeyContext,R::Ring)
    M = free_module(R,1)
    return constant_mackey_functor(context,M)
end