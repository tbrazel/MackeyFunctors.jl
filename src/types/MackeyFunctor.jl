abstract type AbstractShiftedMackeyFunctor end

const ImmutableGeneratorWord = Tuple{Vararg{Tuple{GeneratorIndex,Int}}}

"""
    MackeyFunctor(G)

todo
"""
struct MackeyFunctor
    context::MackeyContext
    values::Vector{AbstractAlgebra.FPModule}
    cover_restrictions::Vector{Generic.ModuleHomomorphism}
    cover_transfers::Vector{Generic.ModuleHomomorphism}
    # generator_conjugations[i ,j] is the conjugation map c_{g_i} : M(H_j) ->  M(g_i H_j g_i^{-1})
    generator_conjugations::Matrix{Generic.ModuleIsomorphism}
    restriction_cache::Dict{Tuple{SubgroupIndex,SubgroupIndex},Generic.ModuleHomomorphism}
    transfer_cache::Dict{Tuple{SubgroupIndex,SubgroupIndex},Generic.ModuleHomomorphism}
    conjugation_cache::Dict{
        Tuple{SubgroupIndex,ImmutableGeneratorWord},
        Generic.ModuleIsomorphism,
    }
    shift_cache::Dict{SubgroupIndex,AbstractShiftedMackeyFunctor}

    function MackeyFunctor(
        context::MackeyContext,
        values::AbstractVector{<:AbstractAlgebra.FPModule},
        cover_restrictions::AbstractVector{<:Generic.ModuleHomomorphism},
        cover_transfers::AbstractVector{<:Generic.ModuleHomomorphism},
        generator_conjugations::AbstractMatrix{<:Generic.ModuleIsomorphism};
        verify::Bool=true
    )
        result = new(
            context,
            values,
            cover_restrictions,
            cover_transfers,
            generator_conjugations,
            Dict{Tuple{SubgroupIndex,SubgroupIndex},Generic.ModuleHomomorphism}(),
            Dict{Tuple{SubgroupIndex,SubgroupIndex},Generic.ModuleHomomorphism}(),
            Dict{
                Tuple{SubgroupIndex,ImmutableGeneratorWord},
                Generic.ModuleIsomorphism,
            }(),
            Dict{
                SubgroupIndex,
                AbstractShiftedMackeyFunctor
            }()
        )

        if verify
            _verify_mackey_functor(result)
        end

        return result
    end
end

function _verify_mackey_functor(mf::MackeyFunctor)
    _verify_values(mf.context, mf.values)
    _verify_cover_maps(
        mf.context,
        mf.values,
        mf.cover_restrictions,
        mf.cover_transfers,
    )
    _verify_generator_conjugations(
        mf.context,
        mf.values,
        mf.generator_conjugations,
    )
    _verify_conjugation_relations(mf)
    _verify_conjugation_cover_compatibility(mf)
    _verify_subgroup_path_independence!(mf)
    _verify_mackey_double_cosets(mf)
    return nothing
end

function _verify_values(
    context::MackeyContext,
    values::AbstractVector{<:AbstractAlgebra.FPModule},
)
    length(values) == length(context.subgroups) ||
        throw(ArgumentError("There must be one value for each subgroup."))

    # A Mackey functor is valued in modules over one fixed coefficient ring.
    # Several later operations assume this without rechecking it:
    # coefficient_ring reads the first value, zero_homomorphism builds matrices
    # over the domain's base ring, and map comparisons evaluate generators
    # inside common module categories.  Check this before the map-domain checks
    # below, so mixed-ring input fails with the actual structural problem.
    coefficient_ring = base_ring(values[begin])
    for i in eachindex(values)
        base_ring(values[i]) == coefficient_ring ||
            throw(ArgumentError("All values of a Mackey functor must be modules over the same base ring. Value $i has base ring $(base_ring(values[i])), but value $(firstindex(values)) has base ring $coefficient_ring."))
    end

    return nothing
end

function _verify_cover_maps(
    context::MackeyContext,
    values::AbstractVector{<:AbstractAlgebra.FPModule},
    cover_restrictions::AbstractVector{<:Generic.ModuleHomomorphism},
    cover_transfers::AbstractVector{<:Generic.ModuleHomomorphism},
)
    covers = context.covers

    length(cover_restrictions) == length(covers) ||
        throw(ArgumentError("There must be one restriction for each cover."))
    length(cover_transfers) == length(covers) ||
        throw(ArgumentError("There must be one transfer for each cover."))

    for (cover_index, (i, j)) in enumerate(covers)
        restriction_map = cover_restrictions[cover_index]
        domain(restriction_map) === values[j] ||
            throw(ArgumentError("Restriction for cover $((i, j)) has the wrong domain."))
        codomain(restriction_map) === values[i] ||
            throw(ArgumentError("Restriction for cover $((i, j)) has the wrong codomain."))

        transfer_map = cover_transfers[cover_index]
        domain(transfer_map) === values[i] ||
            throw(ArgumentError("Transfer for cover $((i, j)) has the wrong domain."))
        codomain(transfer_map) === values[j] ||
            throw(ArgumentError("Transfer for cover $((i, j)) has the wrong codomain."))
    end

    return nothing
end

function _verify_generator_conjugations(
    context::MackeyContext,
    values::AbstractVector{<:AbstractAlgebra.FPModule},
    generator_conjugations::AbstractMatrix{<:Generic.ModuleIsomorphism},
)
    subgroups = context.subgroups
    generators = context.generators
    generator_left_conjugation_matrix = context.generator_left_conjugation_matrix

    size(generator_conjugations) == (length(generators), length(subgroups)) ||
        throw(ArgumentError("There must be one conjugation for each generator and subgroup."))

    for i in eachindex(subgroups), n in eachindex(generators)
        conjugation_map = generator_conjugations[n, i]
        target_index = generator_left_conjugation_matrix[n, i]

        domain(conjugation_map) === values[i] ||
            throw(ArgumentError("Conjugation for generator $n and subgroup $i has the wrong domain."))
        codomain(conjugation_map) === values[target_index] ||
            throw(ArgumentError("Conjugation for generator $n and subgroup $i has the wrong codomain."))
    end

    return nothing
end

function _verify_conjugation_relations(mf::MackeyFunctor)
    context = mf.context

    # If h lies in H, conjugation by h acts trivially on the orbit G/H.  The
    # Mackey functor data therefore has to make c_h the identity on M(H).  It
    # is enough to check this on a generating set for each subgroup H.
    for (i, H) in enumerate(context.subgroups)
        for h in GAP.Globals.GeneratorsOfGroup(H)
            conj_h_H = conjugation(mf, h, i)
            is_identity_module_homomorphism(conj_h_H) ||
                throw(ArgumentError("Conjugation by $h at level $H is not the identity"))
        end
    end

    # The supplied generator conjugations must respect the relators in the
    # chosen presentation of G.  The context stores those relators as words in
    # exactly the same generators used to index generator_conjugations.
    for i in eachindex(context.subgroups), relation_word in context.generator_relations
        conj_by_relation_word = conjugation(mf, relation_word, i)
        is_identity_module_homomorphism(conj_by_relation_word) ||
            throw(ArgumentError("Specified conjugations do not form a valid group action."))
    end

    return nothing
end

function _verify_conjugation_cover_compatibility(mf::MackeyFunctor)
    context = mf.context
    generator_left_conjugation_matrix = context.generator_left_conjugation_matrix

    # Webb's compatibility axioms say that conjugating a cover H < K by a
    # generator g commutes with both the restriction and transfer attached to
    # that cover.  Since conjugation preserves covers, gHg^-1 < gKg^-1 is also
    # one of the context covers; compare the two routes around each square.
    for (cover_index, (i, j)) in enumerate(context.covers)
        for n in eachindex(context.generators)
            gHginvs_index = generator_left_conjugation_matrix[n, i]
            gKginvs_index = generator_left_conjugation_matrix[n, j]
            index_of_conjugated_cover =
                first(context.paths[(gHginvs_index, gKginvs_index)])

            mf.cover_restrictions[cover_index] * mf.generator_conjugations[n, i] ==
            mf.generator_conjugations[n, j] * mf.cover_restrictions[index_of_conjugated_cover] ||
                throw(ArgumentError("Cover restrictions don't commute with generator conjugation."))

            mf.cover_transfers[cover_index] * mf.generator_conjugations[n, j] ==
            mf.generator_conjugations[n, i] * mf.cover_transfers[index_of_conjugated_cover] ||
                throw(ArgumentError("Cover transfers don't commute with generator conjugation."))
        end
    end

    return nothing
end

function _verify_subgroup_path_independence!(mf::MackeyFunctor)
    context = mf.context

    # Build restrictions and transfers along arbitrary subgroup inclusions by
    # induction on chains of covers in the subgroup lattice.
    #
    # The dictionary has one entry for each inclusion H <= K whose composite
    # maps have already been constructed.  Its value is a pair (tr, res), where
    # tr: M(H) -> M(K) and res: M(K) -> M(H).
    #
    # We begin with cover inclusions, since those maps are part of the input.
    # When an inclusion H <= K is taken from the queue, we extend it only across
    # covers K < L.  If H <= L has already appeared, then we have found a second
    # chain of covers from H to L, and we check that the two resulting transfer
    # maps agree and that the two resulting restriction maps agree.
    dictionary_of_paths = Dict{
        Tuple{SubgroupIndex,SubgroupIndex},
        Tuple{Generic.ModuleHomomorphism,Generic.ModuleHomomorphism},
    }()

    for (n, cover) in enumerate(context.covers)
        dictionary_of_paths[cover] = (
            mf.cover_transfers[n],
            mf.cover_restrictions[n],
        )
    end

    outgoing_covers = [Tuple{Int,SubgroupIndex}[] for _ in eachindex(context.subgroups)]
    for (n, (i, j)) in enumerate(context.covers)
        push!(outgoing_covers[i], (n, j))
    end

    queue = collect(keys(dictionary_of_paths))
    queue_head = 1
    while queue_head <= length(queue)
        H_index, K_index = queue[queue_head]
        queue_head += 1

        tr, res = dictionary_of_paths[(H_index, K_index)]
        for (n, next_K_index) in outgoing_covers[K_index]
            new_key = (H_index, next_K_index)
            candidate_tr = tr * mf.cover_transfers[n]
            candidate_res = mf.cover_restrictions[n] * res

            if haskey(dictionary_of_paths, new_key)
                existing_tr, existing_res = dictionary_of_paths[new_key]
                existing_tr == candidate_tr ||
                    throw(ArgumentError("Transfers do not agree along all possible subgroup paths."))
                existing_res == candidate_res ||
                    throw(ArgumentError("Restrictions do not agree along all possible subgroup paths."))
            else
                dictionary_of_paths[new_key] = (candidate_tr, candidate_res)
                push!(queue, new_key)
            end
        end
    end

    for (key, (tr, res)) in dictionary_of_paths
        mf.transfer_cache[key] = tr
        mf.restriction_cache[key] = res
    end

    return nothing
end

function _verify_mackey_double_cosets(mf::MackeyFunctor)
    context = mf.context

    # It suffices to check the double-coset formula for spans J < H > K where
    # both inclusions are covers.  The context has already precomputed the
    # needed representatives and intersections for exactly those triples.
    for (n1, (j, h)) in enumerate(context.covers)
        for (n2, (k, l)) in enumerate(context.covers)
            h == l || continue

            dc_infos = double_coset_infos(context, j, h, k)
            dc_lhs = mf.cover_transfers[n2] * mf.cover_restrictions[n1]
            dc_rhs = zero_homomorphism(domain(dc_lhs), codomain(dc_lhs))

            for info in dc_infos
                JxcapK_index = info.left_conjugate_intersection_index
                JcapxK_index =
                    info.left_intersection_conjugated_right_index

                dc_restriction = restriction(mf, JxcapK_index, k)
                dc_transfer = transfer(mf, JcapxK_index, j)
                dc_conjugation = conjugation(
                    mf,
                    info.representative,
                    JxcapK_index,
                )

                dc_rhs += dc_restriction * dc_conjugation * dc_transfer
            end

            dc_lhs == dc_rhs || throw(ArgumentError("Double coset formula failed."))
        end
    end

    return nothing
end

"""
    value(M,i)

Returns the value of the Mackey functor `M` at the subgroup index `i`.
"""
function value(mf::MackeyFunctor, H_idx::SubgroupIndex)
    return mf.values[H_idx]
end


"""
    restriction(M,i,j)

Given a Mackey functor `M` and two `SubgroupIndex` values `i` and `j` corresponding to subgroups ``H = H[i] \\le H[j] = K \\le G``, this method returns the restriction homomorphism ``M(K) \\to M(H)`` as the type `AbstractAlgebra.Generic.ModuleHomomorphism`.
"""
function restriction(mf::MackeyFunctor, H_index::SubgroupIndex, K_index::SubgroupIndex)
    key = (H_index, K_index)
    if haskey(mf.restriction_cache, key)
        return mf.restriction_cache[key]
    end

    # Make sure H<K first
    is_subgroup(mf.context, H_index, K_index) || throw(ArgumentError("There must exist a path from subgroup 1 to subgroup 2 in order to restrict."))

    path_indices = mf.context.paths[(H_index, K_index)]
    # Start with identity on M(H)
    result = identity_homomorphism(value(mf, H_index))

    #
    for idx in path_indices
        result = mf.cover_restrictions[idx] * result
    end

    mf.restriction_cache[key] = result
    return result
end

"""
    transfer(M,i,j)

Given a Mackey functor `M` and two `SubgroupIndex` values `i` and `j` corresponding to subgroups ``H = H[i] \\le H[j] = K \\le G``, this method returns the transfer homomorphism ``M(H) \\to M(K)`` as the type `AbstractAlgebra.Generic.ModuleHomomorphism`.
"""
function transfer(mf::MackeyFunctor, H_index::SubgroupIndex, K_index::SubgroupIndex)
    key = (H_index, K_index)
    if haskey(mf.transfer_cache, key)
        return mf.transfer_cache[key]
    end

    is_subgroup(mf.context, H_index, K_index) || throw(ArgumentError("There must exist a path from subgroup 1 to subgroup 2 in order to transfer."))
    path_indices = mf.context.paths[(H_index, K_index)]

    result = identity_homomorphism(value(mf, H_index))

    for idx in path_indices
        result = result * mf.cover_transfers[idx]
    end

    mf.transfer_cache[key] = result
    return result
end

"""
    conjugation(M,g,i)

Given a Mackey functor `M`, a `GroupElement` ``g\\in G``, and a `SubgroupIndex` `i`, this method returns the conjugation isomorphism ``M(H[i])\\to M(gH[i]g^{-1})`` as an `AbstractAlgebra.Generic.ModuleIsomorphism`.
"""
function conjugation(mf::MackeyFunctor, g::GroupElement, H_idx::SubgroupIndex)::Generic.ModuleIsomorphism
    conjugation(mf, generator_word(mf.context, g), H_idx)
end

function conjugation(mf::MackeyFunctor, word::GeneratorWord, H_idx::SubgroupIndex)::Generic.ModuleIsomorphism
    key = (H_idx, Tuple(word))
    if haskey(mf.conjugation_cache, key)
        return mf.conjugation_cache[key]
    end

    result = identity_isomorphism(value(mf, H_idx))
    target_of_result = H_idx
    for (g, n) in reverse(word)
        if n>0
            for j in 1:n
                result = result * mf.generator_conjugations[g, target_of_result]

                target_of_result = mf.context.generator_left_conjugation_matrix[g, target_of_result]
            end

        else
            for j in 1:abs(n)
                target_of_result = mf.context.generator_right_conjugation_matrix[g, target_of_result]

                result = result * inv(mf.generator_conjugations[g, target_of_result])
            end
        end
    end
    mf.conjugation_cache[key] = result
    result
end

"""
    coefficient_ring(M)

Returns the underlying coefficient ring of the Mackey functor `M`.
"""
function coefficient_ring(mf::MackeyFunctor)
    return base_ring(mf.values[1])
end

function shift(mf::MackeyFunctor,
    H_index::SubgroupIndex;
    verify::Bool=true,
)::MackeyFunctor
    if haskey(mf.shift_cache, H_index)
        return mf.shift_cache[H_index].underlying_mackey_functor
    end
    # TODO: shifting by G/G should return mf again, but as a *shifted Mackey functor* -- just a speed improvement
    mf.shift_cache[H_index] = _shift(mf, H_index; verify)
    return mf.shift_cache[H_index].underlying_mackey_functor
end
