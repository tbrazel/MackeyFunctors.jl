"""
    direct_sum(M::MackeyFunctor, N::MackeyFunctor) -> (MackeyFunctor, Vector, Vector)

Direct sum of two Mackey functors defined over the same Mackey context.
Returns `(M_plus_N, injections, projections)`, where `injections[1]` is the
canonical map `M -> M_plus_N`, `injections[2]` is `N -> M_plus_N`, and the
projections go in the opposite directions.
"""
function direct_sum(
    M::MackeyFunctor,
    N::MackeyFunctor
)

    M.context === N.context || throw(ArgumentError("Mackey functors are not defined over the same context."))

    context = M.context

    direct_sum_data = [
        direct_sum(AbstractAlgebra.FPModule[
            M.values[subgp_idx],
            N.values[subgp_idx],
        ])
        for subgp_idx in eachindex(M.values)
    ]
    values = AbstractAlgebra.FPModule[first(data) for data in direct_sum_data]

    cover_restrictions = similar(M.cover_restrictions)
    cover_transfers = similar(M.cover_transfers)
    for (cover_index, (i, j)) in enumerate(context.covers)
        p = M.cover_restrictions[cover_index]
        q = N.cover_restrictions[cover_index]
        cover_restrictions[cover_index] = direct_sum_homomorphism(values[j], values[i], p, q)

        p = M.cover_transfers[cover_index]
        q = N.cover_transfers[cover_index]
        cover_transfers[cover_index] = direct_sum_homomorphism(values[i], values[j], p, q)
    end

    generator_conjugations = similar(M.generator_conjugations)
    for i in eachindex(M.context.subgroups), n in eachindex(M.context.generators)
        p = M.generator_conjugations[n, i]
        q = N.generator_conjugations[n, i]
        source = values[i]
        target = values[context.generator_left_conjugation_matrix[n, i]]
        conjugation_hom = direct_sum_homomorphism(source, target, p, q)

        generator_conjugations[n, i] = ModuleIsomorphism(
            source,
            target,
            matrix(conjugation_hom),
        )
    end

    sum_functor = MackeyFunctor(
        context,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations,
        #verify=false,
    )

    left_injection = MackeyFunctorHomomorphism(
        M,
        sum_functor,
        Generic.ModuleHomomorphism[data[2][1] for data in direct_sum_data],
    )
    right_injection = MackeyFunctorHomomorphism(
        N,
        sum_functor,
        Generic.ModuleHomomorphism[data[2][2] for data in direct_sum_data],
    )
    left_projection = MackeyFunctorHomomorphism(
        sum_functor,
        M,
        Generic.ModuleHomomorphism[data[3][1] for data in direct_sum_data],
    )
    right_projection = MackeyFunctorHomomorphism(
        sum_functor,
        N,
        Generic.ModuleHomomorphism[data[3][2] for data in direct_sum_data],
    )

    return (
        sum_functor,
        MackeyFunctorHomomorphism[left_injection, right_injection],
        MackeyFunctorHomomorphism[left_projection, right_projection],
    )
end
