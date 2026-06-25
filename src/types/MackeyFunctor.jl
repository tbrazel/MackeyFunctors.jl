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
        )

        # TODO: sanity check that all the values are modules over the same base ring

        # If the user doesn't want verification, we just return the result. Default setting is to verify
        if !verify
            return result
        end

        # Pull stuff from the Mackey context
        G = context.group
        subgroups = context.subgroups
        covers = context.covers
        generators = context.generators
        generator_left_conjugation_matrix = context.generator_left_conjugation_matrix
        double_coset_formulae = context.double_coset_formulae


        # Before the Mackey axioms, we check that the input is well-formed, i.e. that the values, restrictions, transfers, and conjugations have the right domains and codomains (and that there is the right number of each)
        length(values) == length(subgroups) || throw(ArgumentError("There must be one value for each subgroup."))
        length(cover_restrictions) == length(covers) || throw(ArgumentError("There must be one restriction for each cover."))
        length(cover_transfers) == length(covers) || throw(ArgumentError("There must be one transfer for each cover."))
        size(generator_conjugations) == (length(generators), length(subgroups)) || throw(ArgumentError("There must be one conjugation for each generator and subgroup."))

        for (cover_index, (i, j)) in enumerate(covers)
            restriction_map = cover_restrictions[cover_index]
            domain(restriction_map) === values[j] || throw(ArgumentError("Restriction for cover $((i, j)) has the wrong domain."))
            codomain(restriction_map) === values[i] || throw(ArgumentError("Restriction for cover $((i, j)) has the wrong codomain."))

            transfer_map = cover_transfers[cover_index]
            domain(transfer_map) === values[i] || throw(ArgumentError("Transfer for cover $((i, j)) has the wrong domain."))
            codomain(transfer_map) === values[j] || throw(ArgumentError("Transfer for cover $((i, j)) has the wrong codomain."))
        end

        for i in eachindex(subgroups), n in eachindex(generators)
            conjugation_map = generator_conjugations[n, i]
            target_index = generator_left_conjugation_matrix[n, i]

            domain(conjugation_map) === values[i] || throw(ArgumentError("Conjugation for generator $n and subgroup $i has the wrong domain."))
            codomain(conjugation_map) === values[target_index] || throw(ArgumentError("Conjugation for generator $n and subgroup $i has the wrong codomain."))
        end

        # 1. For h in H, conjugation by h is the identity on M(H)
        for (i, H) in enumerate(subgroups)

            # It suffices to check on generators of H
            for h in GAP.Globals.GeneratorsOfGroup(H)
                conj_h_H = conjugation(result, h, i)
                is_identity_module_homomorphism(conj_h_H) || throw(ArgumentError("Conjugation by $h at level $H is not the identity"))
            end
        end

        # 2. The relations between the generators of G are satisfied by the conjugation automorphisms.
        relation_words = generator_relations_from_isomorphism(context.fp_isomorphism)
        for i in eachindex(subgroups), relation_word in relation_words
            # For each subgroup H, and for every relation (viewed as a word in the generators)
            # We build the map which conjugates M(H) by the relation word
            conj_by_relation_word = conjugation(result, relation_word, i)

            # We assert this is the identity homomorphism on M(H)
            is_identity_module_homomorphism(conj_by_relation_word) || throw(ArgumentError("Specified conjugations do not form a valid group action."))
        end

        # 3. Check Webb Axiom 4 and 5 for covers (compatibility of transfers/conjugation and restriction/conjugation)
        for (cover_index, (i, j)) in enumerate(covers)
            for (n, g) in enumerate(generators)

                # We have res/tr between M(H) and M(K), and we have a generator g. Now we want to check that conjugation by g commutes with res/tr on M(gHg^{-1}) and M(gKg^{-1})

                # First we get the indices of the subgroups gHg^{-1} and gKg^{-1} in our subgroup list
                gHginvs_index = generator_left_conjugation_matrix[n, i]
                gKginvs_index = generator_left_conjugation_matrix[n, j]

                # Since H<K was a cover, we know gHg^{-1}<gKg^{-1} must be a cover, so we get its index
                index_of_conjugated_cover = first(context.paths[(gHginvs_index, gKginvs_index)])

                # We assert restriction along covers commutes with conjugation by generators
                map_eq(
                    cover_restrictions[cover_index] * generator_conjugations[n, i],
                    generator_conjugations[n, j] * cover_restrictions[index_of_conjugated_cover]
                ) || throw(ArgumentError("Cover restrictions don't commute with generator conjugation."))

                # We assert transfer along covers commutes with conjugation by generators
                map_eq(
                    cover_transfers[cover_index] * generator_conjugations[n, j],
                    generator_conjugations[n, i] * cover_transfers[index_of_conjugated_cover]
                ) || throw(ArgumentError("Cover transfers don't commute with generator conjugation."))
            end
        end

        # 4. For H<K not a cover, check any composite of covers beginning at H and ending at K yields the same well-defined transfer and restriction

        # Strategy: we build the restrictions and transfers along arbitrary subgroup
        # inclusions by induction on chains of covers in the subgroup lattice.
        #
        # The dictionary has one entry for each inclusion H <= K whose composite maps
        # we have already constructed. Its value is a pair (tr, res), where
        # tr: M(H) -> M(K) and res: M(K) -> M(H).
        #
        # We begin with the cover inclusions, since those maps are part of the input.
        # Then we keep a queue of the inclusions whose composites have just become
        # known. When an inclusion H <= K is taken from the queue, we extend it only
        # across covers K < L. This gives a candidate composite for H <= L.
        #
        # If H <= L has not appeared before, we record this candidate and add H <= L
        # to the queue, so it can itself be extended later. If H <= L has appeared
        # before, then we have found a second chain of covers from H to L, and we
        # check that the two resulting transfer maps agree and that the two resulting
        # restriction maps agree.
        #
        # Since the subgroup lattice is finite, this process eventually considers
        # every composite of cover maps. Thus the input cover maps determine
        # well-defined restrictions and transfers along all subgroup inclusions
        # exactly when every repeated inclusion gives the same maps.
        dictionary_of_paths = Dict{Tuple{SubgroupIndex,SubgroupIndex},Tuple{Generic.ModuleHomomorphism,Generic.ModuleHomomorphism}}()

        # We first initialize the dictionaries with restriction and transfer along covers
        for (n, cov) in enumerate(context.covers)
            dictionary_of_paths[cov] = (cover_transfers[n], cover_restrictions[n])
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
                # Extend the path H < K by a cover K < next_K.
                # tr goes M(H) -> M(K)
                # res goes M(K) -> M(H)

                new_key = (H_index, next_K_index)

                candidate_tr = tr * cover_transfers[n]
                candidate_res = cover_restrictions[n] * res

                if haskey(dictionary_of_paths, new_key)
                    existing_tr, existing_res = dictionary_of_paths[new_key]
                    map_eq(
                        existing_tr, candidate_tr
                    ) || throw(ArgumentError("Transfers do not agree along all possible subgroup paths."))

                    map_eq(
                        existing_res, candidate_res
                    ) || throw(ArgumentError("Restrictions do not agree along all possible subgroup paths."))

                else
                    dictionary_of_paths[new_key] = (candidate_tr, candidate_res)
                    push!(queue, new_key)
                end
            end
        end

        for (key, (tr, res)) in dictionary_of_paths
            result.transfer_cache[key] = tr
            result.restriction_cache[key] = res
        end

        # 5. Check double coset formula for covers
        for (n1, (j, h)) in enumerate(covers)
            for (n2, (k, l)) in enumerate(covers)
                h == l || continue

                # We have J<H and K<H, which are both covers
                # J = subgroups[j]
                # H = subgroups[h]
                # K = subgroups[k]

                # Pull the double coset representatives for J\H/K
                dc_reps = double_coset_formulae[(j, h, k)]

                # The LHS of the double coset formula is res_J^H tr_K^H
                dc_lhs = cover_transfers[n2] * cover_restrictions[n1]

                # We start with the RHS being the zero map from M(K) -> M(J)
                dc_rhs = zero_homomorphism(domain(dc_lhs), codomain(dc_lhs))

                # For each pair of (g,J^x \cap K), we add the needed term to the RHS of the double coset formula
                for (w, JxcapK_index) in dc_reps
                    JcapxK_index = conjugate_subgroup_by_word(context, JxcapK_index, w)

                    dc_restriction = restriction(result, JxcapK_index, k)
                    dc_transfer = transfer(result, JcapxK_index, j)
                    dc_conjugation = conjugation(result, w, JxcapK_index)

                    dc_rhs += dc_restriction * dc_conjugation * dc_transfer
                end

                # Assert that the LHS and RHS of the double coset formula agree
                map_eq(
                    dc_lhs,
                    dc_rhs
                ) || throw(ArgumentError("Double coset formula failed."))
            end
        end

        return result
    end
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

function Base.show(io::IO, obj::MackeyFunctor)
    println(io, "MackeyFunctor for group ", String(GAP.Globals.StructureDescription(obj.context.group)), " over base ring ", coefficient_ring(obj))
    
    for (i,(h,k)) in enumerate(obj.context.covers)
        hname = String(GAP.Globals.StructureDescription(obj.context.subgroups[h]))
        kname = String(GAP.Globals.StructureDescription(obj.context.subgroups[k]))
        println("Restriction ", kname, " ---> ", hname)
        display(matrix(obj.cover_restrictions[i]))
        println("Transfer ", hname, " ---> ", kname)
        display(matrix(obj.cover_transfers[i]))
    end

end
