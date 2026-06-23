using AbstractAlgebra

struct MackeyFunctor
    context::MackeyContext
    values::Vector{AbstractAlgebra.FPModule}
    cover_restrictions::Vector{Generic.ModuleHomomorphism}
    cover_transfers::Vector{Generic.ModuleHomomorphism}
    # generator_conjugations[i,j] is the conjugation map c_{g_i} : M(H_j) ->  M(g_i H_j g_i^{-1})
    generator_conjugations::Matrix{Generic.ModuleIsomorphism}

    function MackeyFunctor(
        context::MackeyContext,
        values::Vector{AbstractAlgebra.FPModule},
        cover_restrictions::Vector{Generic.ModuleHomomorphism},
        cover_transfers::Vector{Generic.ModuleHomomorphism},
        generator_conjugations::Matrix{Generic.ModuleIsomorphism},
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
                conj_h_H = conjugation(result, h, i)
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
        for (cover_index,(i,j)) in enumerate(covers)
            for (n,g) in enumerate(generators)
                H = subgroups[i]
                K = subgroups[j]

                gHginvs_index = generatorLeftConjugationMatrix[n,i]

                gKginvs_index = generatorLeftConjugationMatrix[n,j]
                
                index_of_conjugated_cover = first(context.paths[(gHginvs_index,gKginvs_index)])
                
                # Check res commutes
                same_module_map(
                    generator_conjugations[n,i]*cover_restrictions[cover_index],
                    cover_restrictions[index_of_conjugated_cover]*generator_conjugations[n,j]
                    ) || throw(ArgumentError("Cover restrictions don't commute with generator conjugation."))

                # Check tr commutes
                same_module_map(
                    generator_conjugations[n,j]*cover_transfers[cover_index],
                    cover_transfers[index_of_conjugated_cover]*generator_conjugations[n,i]
                ) || throw(ArgumentError("Cover transfers don't commute with generator conjugation."))
            end
        end

        # Next, we need to check the restrictions and transfers. This entails:
        # 4. Check double coset formula for covers
        for (n1,(j,h)) in enumerate(covers)
            for (n2,(k,l)) in enumerate(covers)
                h == l || continue

                # We have J<H and K<H both covers
                J = subgroups[j]
                H = subgroups[h]
                K = subgroups[k]
                
                dc_reps = doubleCosetRepresentatives[(j,h,k)]

                dc_lhs = cover_restrictions[n1]*cover_transfers[n2]

                dc_rhs = zero_homomorphism(domain(dc_lhs), codomain(dc_lhs))

                for (w,JxcapK_index) in dc_reps
                    JcapxK_index = conjugate_subgroup_by_word(context,JxcapK_index,w)

                    dc_restriction = restriction(result,JxcapK_index,k)
                    dc_transfer = transfer(result,JcapxK_index, j)
                    dc_conjugation = conjugation(result,JxcapK_index,w)

                    dc_rhs += dc_transfer*dc_conjugation*dc_restriction
                end
                same_module_map(
                    dc_lhs,
                    dc_rhs
                ) || throw(ArgumentError("Double coset formula failed."))
            end
        end

        # 5. For H<K not a cover, check any composite of covers beginning at H and ending at K yields the same well-defined transfer and restriction
        

        return result
    end
end

# Returns the value of a Mackey functor at subgroup index H_idx
function value(mf::MackeyFunctor, H_idx::SubgroupIndex)
    return mf.values[H_idx]
end

# Get the value of restriction M(K) -> M(H) for an arbitrary subgroup inclusion H<K
function restriction(mf::MackeyFunctor,H_index::SubgroupIndex, K_index::SubgroupIndex)
    path_indices = mf.context.paths[(H_index,K_index)]
    # Start with identity on M(H)
    result = identity_isomorphism(value(mf,H_index))

    # 
    for idx in path_indices
        result = result * mf.cover_restrictions[idx] 
    end

    return result
end

function transfer(mf::MackeyFunctor, H_index::SubgroupIndex, K_index::SubgroupInde)
    path_indices = mf.context.paths[(H_index,K_index)]

    result = identity_isomorphism(value(mf,H_index))

    for idx in path.indices
        result = mf.cover_transfers[idx] * result
    end

    return result
end

# IN PROGRESS
function conjugation(mf::MackeyFunctor, H_idx::SubgroupIndex, g::GroupElement)::Generic.ModuleIsomorphism
    G = mf.context.group
    result = identity_isomorphism(value(mf, H_idx))
    target_of_result = H_idx
    word = generator_word(G, g)

    # TODO come back and redo this with conjugate_subgroup_by_word method ?
    for (g, n) in reverse(word)
        if n>0
            for j in 1:n
                result = mf.generator_conjugations[g, target_of_result] * result

                target_of_result = mf.context.generatorLeftConjugationMatrix[g, target_of_result]
            end

        else
            for j in 1::abs(n)
                target_of_result = mf.context.generatorRightConjugationMatrix[g, target_of_result]

                result = inv(mf.generator_conjugations[g, target_of_result]) * result
            end
        end
    end
    result
end

function conjugation(mf::MackeyFunctor, H_idx::SubgroupIndex, word::GeneratorWord)::Generic.ModuleIsomorphism
    G = mf.context.group
    result = identity_isomorphism(value(mf, H_idx))
    target_of_result = H_idx
    # TODO come back and redo this with conjugate_subgroup_by_word method ?
    for (g, n) in reverse(word)
        if n>0
            for j in 1:n
                result = mf.generator_conjugations[g, target_of_result] * result

                target_of_result = mf.context.generatorLeftConjugationMatrix[g, target_of_result]
            end

        else
            for j in 1::abs(n)
                target_of_result = mf.context.generatorRightConjugationMatrix[g, target_of_result]

                result = inv(mf.generator_conjugations[g, target_of_result]) * result
            end
        end
    end
    result
end