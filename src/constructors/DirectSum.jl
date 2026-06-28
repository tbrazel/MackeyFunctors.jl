"""
    direct_sum(mfv::AbstractVector{MackeyFunctor}) -> (MackeyFunctor, Vector, Vector)

Direct sum of the given Mackey functors defined over the same Mackey context.
Returns `(mf_sum, injections, projections)`, where `injections[i]` is the
canonical injection of the `i`-th factor and `projections[i]` the `i`-th projection.
"""
function direct_sum(mfv::AbstractVector{MackeyFunctor})
    isempty(mfv) && throw(ArgumentError("Direct sum of no Mackey functors is not supported."))
    context = mfv[1].context
    all(mf -> mf.context === context, mfv) || throw(ArgumentError("Mackey functors are not defined over the same context."))

    direct_sum_data = [direct_sum(map(mf -> mf.values[i], mfv)) for i in eachindex(context.subgroups)]
    values = map(first, direct_sum_data)

    cover_restrictions = similar(mfv[1].cover_restrictions)
    cover_transfers = similar(mfv[1].cover_transfers)
    for (cover_index, (i, j)) in enumerate(context.covers)
        summand_restrictions = map(mf -> mf.cover_restrictions[cover_index], mfv)
        cover_restrictions[cover_index] = direct_sum(summand_restrictions; domain = values[j], codomain = values[i])

        summand_transfers = map(mf -> mf.cover_transfers[cover_index], mfv)
        cover_transfers[cover_index] = direct_sum(summand_transfers; domain = values[i], codomain = values[j])
    end

    generator_conjugations = Generic.ModuleIsomorphism[
        direct_sum(
            map(mf -> mf.generator_conjugations[n, i], mfv);
            domain = values[i],
            codomain = values[context.generator_left_conjugation_matrix[n, i]],
        ) for n in eachindex(context.generators), i in eachindex(context.subgroups)
    ]

    mf_sum = MackeyFunctor(
        context,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations
    )

    injections = [
        MackeyFunctorHomomorphism(mf, mf_sum, Generic.ModuleHomomorphism[data[2][i] for data in direct_sum_data])
            for (i, mf) in enumerate(mfv)
    ]
    projections = [
        MackeyFunctorHomomorphism(mf_sum, mf, Generic.ModuleHomomorphism[data[3][i] for data in direct_sum_data])
            for (i, mf) in enumerate(mfv)
    ]

    return mf_sum, injections, projections
end

direct_sum(mf::MackeyFunctor...) = direct_sum(collect(mf))

function _compose_mackey_functor_homomorphism_pair(
    first::MackeyFunctorHomomorphism,
    second::MackeyFunctorHomomorphism,
)::MackeyFunctorHomomorphism
    first.codomain === second.domain ||
        throw(ArgumentError("Mackey functor homomorphisms cannot be composed."))

    return MackeyFunctorHomomorphism(
        first.domain,
        second.codomain,
        Generic.ModuleHomomorphism[
            first.components[i] * second.components[i]
            for i in eachindex(first.context.subgroups)
        ],
    )
end

function _direct_sum_mackey_functors(
    summands::AbstractVector{MackeyFunctor},
)
    isempty(summands) &&
        throw(ArgumentError("Cannot take the direct sum of no Mackey functors."))

    sum_functor = first(summands)
    injections = MackeyFunctorHomomorphism[id_homomorphism(sum_functor)]
    projections = MackeyFunctorHomomorphism[id_homomorphism(sum_functor)]

    for next_summand in summands[2:end]
        new_sum, binary_injections, binary_projections =
            direct_sum(sum_functor, next_summand)
        left_injection, right_injection = binary_injections
        left_projection, right_projection = binary_projections

        injections = MackeyFunctorHomomorphism[
            _compose_mackey_functor_homomorphism_pair(injection, left_injection)
            for injection in injections
        ]
        push!(injections, right_injection)

        projections = MackeyFunctorHomomorphism[
            _compose_mackey_functor_homomorphism_pair(left_projection, projection)
            for projection in projections
        ]
        push!(projections, right_projection)

        sum_functor = new_sum
    end

    return sum_functor, injections, projections
end

function _check_mackey_functor_homomorphism_matrix(
    maps::AbstractMatrix{<:MackeyFunctorHomomorphism},
)
    isempty(maps) &&
        throw(ArgumentError("Cannot build a homomorphism from an empty matrix."))

    for row in axes(maps, 1), column in axes(maps, 2)
        if maps[row, column].domain !== maps[row, first(axes(maps, 2))].domain
            throw(ArgumentError("Entries in each matrix row must have the same domain."))
        end
        if maps[row, column].codomain !== maps[first(axes(maps, 1)), column].codomain
            throw(ArgumentError("Entries in each matrix column must have the same codomain."))
        end
    end

    return nothing
end

"""
    block_homomorphism(maps::AbstractMatrix{<:MackeyFunctorHomomorphism})

Build the homomorphism represented by a matrix of Mackey functor homomorphisms.
This follows AbstractAlgebra's matrix convention: rows index domain summands
and columns index codomain summands.  Thus `maps[i, j]` must be a map from the
`i`th domain summand to the `j`th codomain summand, and the result has type
`direct_sum(row domains...) -> direct_sum(column codomains...)`.

In particular, a `1 x n` row matrix gives a map from one Mackey functor into a
direct sum of `n` codomains, while an `m x 1` column matrix gives a map from a
direct sum of `m` domains into one Mackey functor.
"""
function block_homomorphism(
    maps::AbstractMatrix{<:MackeyFunctorHomomorphism},
)::MackeyFunctorHomomorphism
    _check_mackey_functor_homomorphism_matrix(maps)

    row_indices = collect(axes(maps, 1))
    column_indices = collect(axes(maps, 2))
    domain_summands = MackeyFunctor[
        maps[row, first(column_indices)].domain
        for row in row_indices
    ]
    codomain_summands = MackeyFunctor[
        maps[first(row_indices), column].codomain
        for column in column_indices
    ]

    domain_sum, _, domain_projections =
        _direct_sum_mackey_functors(domain_summands)
    codomain_sum, codomain_injections, =
        _direct_sum_mackey_functors(codomain_summands)

    components = Generic.ModuleHomomorphism[]
    for H_index in eachindex(domain_sum.context.subgroups)
        component = zero_homomorphism(
            domain_sum.values[H_index],
            codomain_sum.values[H_index],
        )

        for (row_number, row) in enumerate(row_indices), (column_number, column) in enumerate(column_indices)
            component +=
                domain_projections[row_number].components[H_index] *
                maps[row, column].components[H_index] *
                codomain_injections[column_number].components[H_index]
        end

        push!(components, component)
    end

    return MackeyFunctorHomomorphism(domain_sum, codomain_sum, components)
end

function block_homomorphism(
    maps::AbstractVector{<:MackeyFunctorHomomorphism};
    orientation::Symbol=:column,
)::MackeyFunctorHomomorphism
    if orientation === :column
        return block_homomorphism(reshape(collect(maps), :, 1))
    elseif orientation === :row
        return block_homomorphism(reshape(collect(maps), 1, :))
    else
        throw(ArgumentError("orientation must be either :row or :column."))
    end
end

MackeyFunctorHomomorphism(
    maps::AbstractMatrix{<:MackeyFunctorHomomorphism},
) = block_homomorphism(maps)

function MackeyFunctorHomomorphism(
    maps::AbstractVector{<:MackeyFunctorHomomorphism};
    orientation::Symbol=:column,
)
    return block_homomorphism(maps; orientation=orientation)
end

function _zero_mackey_functor_homomorphism(
    domain_mf::MackeyFunctor,
    codomain_mf::MackeyFunctor,
)::MackeyFunctorHomomorphism
    domain_mf.context === codomain_mf.context ||
        throw(ArgumentError("Mackey functors must have the same context."))

    return MackeyFunctorHomomorphism(
        domain_mf,
        codomain_mf,
        Generic.ModuleHomomorphism[
            zero_homomorphism(domain_mf.values[i], codomain_mf.values[i])
            for i in eachindex(domain_mf.context.subgroups)
        ],
    )
end

"""
    direct_sum(f::MackeyFunctorHomomorphism, g::MackeyFunctorHomomorphism)

Return the block-diagonal direct sum homomorphism `f ⊕ g`.
If `f : A -> C` and `g : B -> D`, the result is a map
`A ⊕ B -> C ⊕ D`.
"""
function direct_sum(
    f::MackeyFunctorHomomorphism,
    g::MackeyFunctorHomomorphism,
)::MackeyFunctorHomomorphism
    return block_homomorphism([
        f _zero_mackey_functor_homomorphism(f.domain, g.codomain)
        _zero_mackey_functor_homomorphism(g.domain, f.codomain) g
    ])
end
