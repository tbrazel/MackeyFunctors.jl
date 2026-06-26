"""
    direct_sum(M::MackeyFunctor, N::MackeyFunctor) -> MackeyFunctor

Direct sum of two Mackey functors defined over the same Mackey context.
"""
function direct_sum(
    M::MackeyFunctor,
    N::MackeyFunctor
)

    M.context === N.context || throw(ArgumentError("Mackey functors are not defined over the same context."))

    context = M.context

    values = [
        _direct_sum_module(AbstractAlgebra.FPModule[
            M.values[subgp_idx],
            N.values[subgp_idx],
        ])
        for subgp_idx in eachindex(M.values)
    ]

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

    return MackeyFunctor(
        context,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations,
        #verify=false,
    )
end
