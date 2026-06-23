using MackeyContext
using AbstractAlgebra

struct MackeyFunctor
    context::MackeyContext
    values::Vector{FPModule}
    cover_restrictions::Vector{Generic.ModuleHomomorphism}
    cover_transfers::Vector{Generic.ModuleHomomorphism}
    # generator_conjugations[i,j] is the conjugation map c_{g_i} : M(H_j) ->  M(g_i H_j g_i^{-1})
    generator_conjugations::Matrix{Generic.ModuleIsomorphism}

    function MackeyFunctor(
        context::MackeyContext,
        values::Vector{FpModule},
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
        # 3. Check Webb Axiom 4 and 5 for covers (compatibililty of transfers/conjugation and restriction/conjugation)

        # Next, we need to check the restrictions and transfers. This entails:
        # 4. Check double coset formula for covers

        # 5. For H<K not a cover, check any composite of covers beginning at H and ending at K yields the same well-defined transfer and restriction






        # do stuff here

        return result
    end
end

# Returns the value of a Mackey functor at subgroup index H_idx
function value(mf::MackeyFunctor, H_idx::SubgroupIndex)
    return mf.values[H_idx]
end

# Given a module M, returns its identity as a type Generic.ModuleIsomorphism
function identity_isomorphism(M::AbstractAlgebra.FPModule)::Generic.ModuleIsomorphism
    ModuleIsomorphism(M, M, identity_matrix(base_ring(M), ngens(M)))
end

function is_zero_module_homomorphism(phi::Generic.ModuleHomomorphism)
    return all(x -> phi(x) == 0, gens(domain(phi)))
end

function is_equal_module_homomorphism(phi::Generic.ModuleHomomorphism, psi::Generic.ModuleHomomorphism)
    return domain(phi) === domain(psi) && codomain(phi) === codomain(psi) && is_zero_module_homomorphism(phi - psi)
end

function is_identity_module_homomorphism(phi::Generic.ModuleIsomorphism)
    return domain(phi) === codomain(phi) && all(phi(x) == identity_isomorphism(domain(phi))(x) for x in gens(domain(phi)))
end

function same_module_map(f, g)
    domain(f) === domain(g) || return false
    codomain(f) === codomain(g) || return false
    return all(x -> f(x) == g(x), gens(domain(f)))
end

# IN PROGRESS
function conjugation(mf::MackeyFunctor, H_idx::SubgroupIndex, g::GroupElement)::Generic.ModuleIsomorphism
    G = mf.context.group
    H = mf.context.subgroups[H_idx]
    result = identity_isomorphism(value(mf, H_idx))
    target_of_result = H_idx
    word = generator_word(G, g)


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