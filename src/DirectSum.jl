using AbstractAlgebra

"""
Direct sum of two Mackey functors defined over the same Mackey context.
"""
function direct_sum_mf(
    M::MackeyFunctor,
    N::MackeyFunctor
)

    same_context(M.context, N.context) || throw(ArgumentError("Mackey functors are not defined over the same context."))

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

function direct_sum_homomorphism(
    source,
    target,
    p::Generic.ModuleIsomorphism,
    q::Generic.ModuleIsomorphism,
)
    return direct_sum_homomorphism(source, target, as_homomorphism(p), as_homomorphism(q))
end

function direct_sum_homomorphism(
    source,
    target,
    p::Generic.ModuleHomomorphism,
    q::Generic.ModuleHomomorphism,
)
    R = base_ring(source)
    all(
        mod -> base_ring(mod) == R,
        (target, domain(p), codomain(p), domain(q), codomain(q)),
    ) || throw(ArgumentError("All direct-sum homomorphism modules must have the same base ring."))

    source_generators = ngens(domain(p)) + ngens(domain(q))
    target_generators = ngens(codomain(p)) + ngens(codomain(q))
    ngens(source) == source_generators ||
        throw(ArgumentError("The source does not have the expected direct-sum presentation."))
    ngens(target) == target_generators ||
        throw(ArgumentError("The target does not have the expected direct-sum presentation."))

    block_matrix = zero_matrix(R, ngens(source), ngens(target))
    _copy_matrix_block!(block_matrix, matrix(p), 0, 0)
    _copy_matrix_block!(
        block_matrix,
        matrix(q),
        ngens(domain(p)),
        ngens(codomain(p)),
    )

    return ModuleHomomorphism(source, target, block_matrix)
end

function _copy_matrix_block!(target_matrix, source_matrix, row_offset::Int, column_offset::Int)
    for row in 1:nrows(source_matrix), column in 1:ncols(source_matrix)
        target_matrix[row_offset + row, column_offset + column] =
            source_matrix[row, column]
    end

    return target_matrix
end

function same_context(a::MackeyContext, b::MackeyContext)
    return all(f -> getfield(a, f) == getfield(b, f), fieldnames(MackeyContext))
end
