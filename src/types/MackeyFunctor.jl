# using AbstractAlgebra


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

    function MackeyFunctor(
        context::MackeyContext,
        values::Vector{<:AbstractAlgebra.FPModule},
        cover_restrictions::Vector{<:Generic.ModuleHomomorphism},
        cover_transfers::Vector{<:Generic.ModuleHomomorphism},
        generator_conjugations::Matrix{<:Generic.ModuleIsomorphism},
        verify::Bool=true
    )
        result = new(
            context,
            values,
            cover_restrictions,
            cover_transfers,
            generator_conjugations,
        )

        # If the user doesn't want verification, ok then I guess
        if !verify
            return result
        end

        # Pull stuff from the Mackey context
        G = context.group
        subgroups = context.subgroups
        covers = context.covers
        generators = context.generators
        generatorLeftConjugationMatrix = context.generatorLeftConjugationMatrix
        generatorRightConjugationMatrix = context.generatorRightConjugationMatrix
        doubleCosetRepresentatives = context.doubleCosetRepresentatives


        # First thing to verify: conjugations are valid. This means:
        # 1. For H <= G and h \in H, conjugation by h at level G/H should be identity
        for (i, H) in enumerate(subgroups)
            for h in GAP.Globals.GeneratorsOfGroup(H)
                conj_h_H = conjugation(result, i, h)
                is_identity_module_homomorphism(conj_h_H) || throw(ArgumentError("Conjugation by $h at level $H is not the identity"))
            end
        end
        # 2. The relations between the generators of G are satisfied by the conjugation automorphisms.
        for i in eachindex(subgroups)
            for relation_word in generator_relations(G, generators)
                conj_by_relation_word = conjugation(result, i, relation_word)

                is_identity_module_homomorphism(conj_by_relation_word) || throw(ArgumentError("Specified conjugations do not form a valid group action."))
            end
        end

        # 3. Check Webb Axiom 4 and 5 for covers (compatibililty of transfers/conjugation and restriction/conjugation)
        for (cover_index, (i, j)) in enumerate(covers)
            for (n, g) in enumerate(generators)
                H = subgroups[i]
                K = subgroups[j]

                gHginvs_index = generatorLeftConjugationMatrix[n, i]

                gKginvs_index = generatorLeftConjugationMatrix[n, j]

                index_of_conjugated_cover = first(context.paths[(gHginvs_index, gKginvs_index)])

                # Check res commutes
                map_eq(
                    cover_restrictions[cover_index] * generator_conjugations[n, i],
                    generator_conjugations[n, j] * cover_restrictions[index_of_conjugated_cover]
                ) || throw(ArgumentError("Cover restrictions don't commute with generator conjugation."))

                # Check tr commutes
                map_eq(
                    cover_transfers[cover_index] * generator_conjugations[n, j],
                    generator_conjugations[n, i] * cover_transfers[index_of_conjugated_cover]
                ) || throw(ArgumentError("Cover transfers don't commute with generator conjugation."))
            end
        end

        # Next, we need to check the restrictions and transfers. This entails:
        # 4. Check double coset formula for covers
        for (n1, (j, h)) in enumerate(covers)
            for (n2, (k, l)) in enumerate(covers)
                h == l || continue

                # We have J<H and K<H both covers
                # J = subgroups[j]
                # H = subgroups[h]
                # K = subgroups[k]

                dc_reps = doubleCosetRepresentatives[(j, h, k)]

                dc_lhs = cover_transfers[n2] * cover_restrictions[n1]

                dc_rhs = zero_homomorphism(domain(dc_lhs), codomain(dc_lhs))

                for (w, JxcapK_index) in dc_reps
                    JcapxK_index = conjugate_subgroup_by_word(context, JxcapK_index, w)

                    dc_restriction = restriction(result, JxcapK_index, k)
                    dc_transfer = transfer(result, JcapxK_index, j)
                    dc_conjugation = conjugation(result, JxcapK_index, w)

                    dc_rhs += dc_restriction * dc_conjugation * dc_transfer
                end
                map_eq(
                    dc_lhs,
                    dc_rhs
                ) || throw(ArgumentError("Double coset formula failed."))
            end
        end

        # 5. For H<K not a cover, check any composite of covers beginning at H and ending at K yields the same well-defined transfer and restriction

        # Dictionary has keys (i,j) corresponding to H[i]<H[j], and values (t,r) for t: M(H[i]) -> M(H[j]) a transfer, and r: M(H[j]) -> M(H[i]) a restriction
        dictionary_of_paths = Dict{Tuple{SubgroupIndex,SubgroupIndex},Tuple{Generic.ModuleHomomorphism,Generic.ModuleHomomorphism}}()


        # Initialize the dictionaries with restriction and transfer along covers
        for (n, cov) in enumerate(context.covers)
            dictionary_of_paths[cov] = (cover_transfers[n], cover_restrictions[n])
        end

        # Iterate and see 
        changed = true
        while changed
            changed = false

            for (n, (i, j)) in enumerate(context.covers)
                for ((H_index, K_index), (tr, res)) in dictionary_of_paths
                    K_index == i || continue

                    # So path goes H<K = H[i] < H[j]
                    # tr goes M(H) -> M(K)
                    # res goes M(K) -> M(H)

                    new_key = (H_index, j)

                    candidate_tr = tr * cover_transfers[n]
                    candidate_res = cover_restrictions[n] * res
                    candidate_value = (
                        candidate_tr, candidate_res
                    )

                    if haskey(dictionary_of_paths, new_key)
                        (existing_tr, existing_res) = dictionary_of_paths[new_key]
                        map_eq(
                            existing_tr, candidate_tr
                        ) || throw(ArgumentError("Transfers do not agree along all possible subgroup paths."))

                        map_eq(
                            existing_res, candidate_res
                        ) || throw(ArgumentError("Restrictions do not agree along all possible subgroup paths."))

                    else
                        push!(dictionary_of_paths, new_key=>candidate_value)
                        changed = true
                    end

                end
            end

        end

        return result
    end
end

# Returns the value of a Mackey functor at subgroup index H_idx
function value(mf::MackeyFunctor, H_idx::SubgroupIndex)
    return mf.values[H_idx]
end


"""
    restriction(M,i,j)

If ``H`` is the `i`th subgroup and ``K`` is the `j`th subgroup, this returns the restriction map ``M(K) \\to M(H)``.
"""
function restriction(mf::MackeyFunctor, H_index::SubgroupIndex, K_index::SubgroupIndex)
    # Make sure H<K first
    is_subgroup(mf.context, H_index, K_index) || throw(ArgumentError("There must exist a path from subgroup 1 to subgroup 2 in order to restrict."))

    path_indices = mf.context.paths[(H_index, K_index)]
    # Start with identity on M(H)
    result = identity_homomorphism(value(mf, H_index))

    # 
    for idx in path_indices
        result = mf.cover_restrictions[idx] * result
    end

    return result
end

"""
    transfer(M,H,K)

todo
"""
function transfer(mf::MackeyFunctor, H_index::SubgroupIndex, K_index::SubgroupIndex)
    is_subgroup(mf.context, H_index, K_index) || throw(ArgumentError("There must exist a path from subgroup 1 to subgroup 2 in order to transfer."))
    path_indices = mf.context.paths[(H_index, K_index)]

    result = identity_homomorphism(value(mf, H_index))

    for idx in path_indices
        result = result * mf.cover_transfers[idx]
    end

    return result
end

"""
    conjugation(M,n,g)

If ``H`` denotes the ``n``th subgroup for `G = M.group`, and ``g\\in G`` is a group element, this method returns the conjugation map ``M(H) \\to M(gHg^{-1})`` in the Mackey functor.
"""
function conjugation(mf::MackeyFunctor, H_idx::SubgroupIndex, g::GroupElement)::Generic.ModuleIsomorphism
    conjugation(mf, H_idx, generator_word(mf.context.group, g))
end

function conjugation(mf::MackeyFunctor, H_idx::SubgroupIndex, word::GeneratorWord)::Generic.ModuleIsomorphism
    G = mf.context.group
    result = identity_isomorphism(value(mf, H_idx))
    target_of_result = H_idx
    for (g, n) in reverse(word)
        if n>0
            for j in 1:n
                result = result * mf.generator_conjugations[g, target_of_result]

                target_of_result = mf.context.generatorLeftConjugationMatrix[g, target_of_result]
            end

        else
            for j in 1:abs(n)
                target_of_result = mf.context.generatorRightConjugationMatrix[g, target_of_result]

                result = result * inv(mf.generator_conjugations[g, target_of_result])
            end
        end
    end
    result
end
