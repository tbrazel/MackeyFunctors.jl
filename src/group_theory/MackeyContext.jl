const SubgroupIndex = Int
const GeneratorIndex = Int
const Group = GapObj
const GroupElement = GapObj
const GeneratorWord = Vector{Tuple{GeneratorIndex,Int}}

"""
    MackeyContext(G)

Precompute subgroup-lattice and double-coset data for a finite GAP group `G`.

The field `doubleCosetRepresentatives[(h, k, j)]` is defined for triples where
both `subgroups[h] < subgroups[k]` and `subgroups[j] < subgroups[k]` are cover
relations. Its value is a vector of representatives for
`subgroups[h] \\ subgroups[k] / subgroups[j]`, where each representative is
stored as a word in `generators`.

A word is a vector of `(generator_index, exponent)` pairs. For example,
`[(1, 3), (2, -4)]` represents `generators[1]^3 * generators[2]^-4`, and the
empty vector represents the identity element.
"""
struct MackeyContext
    G::Group
    subgroups::Vector{Group}
    covers::Vector{Tuple{SubgroupIndex,SubgroupIndex}}
    epimorphismFromFreeGroup::GapObj
    generators::Vector{GroupElement}
    generatorLeftConjugationMatrix::Matrix{SubgroupIndex}
    generatorRightConjugationMatrix::Matrix{SubgroupIndex}
    doubleCosetRepresentatives::Dict{
        Tuple{SubgroupIndex,SubgroupIndex,SubgroupIndex},
        Vector{GeneratorWord},
    }

    # Build a MackeyContext from a group G
    function MackeyContext(G::Group)
        GAP.Globals.IsGroup(G) || throw(ArgumentError("Input G must be a GAP group"))

        # get subgroups
        subgroups = Vector{Group}(GAP.Globals.AllSubgroups(G))

        # Build list of covers
        covers = Tuple{SubgroupIndex,SubgroupIndex}[]

        for i in eachindex(subgroups), j in eachindex(subgroups)
            H = subgroups[i]
            K = subgroups[j]

            # GAP.Globals.IsSubgroup(K, H) means H <= K
            H_properly_contained_in_K =
                i != j && Bool(GAP.Globals.IsSubgroup(K, H))

            H_properly_contained_in_K || continue

            has_intermediate = any(eachindex(subgroups)) do l
                L = subgroups[l]

                H_properly_contained_in_L =
                    l != i && Bool(GAP.Globals.IsSubgroup(L, H))

                L_properly_contained_in_K =
                    l != j && Bool(GAP.Globals.IsSubgroup(K, L))

                H_properly_contained_in_L && L_properly_contained_in_K
            end

            has_intermediate || push!(covers, (i, j))
        end
        
        # Get list of generators
        generators = Vector{GroupElement}(GAP.Globals.GeneratorsOfGroup(G))

        # Get epimorphism from free group
        epi_from_free_group = GAP.Globals.EpimorphismFromFreeGroup(G)

        # Build conjugation matrices
        num_rows_conj_matrix = length(subgroups)
        num_cols_conj_matrix = length(generators)
        left_conj_matx =
            Matrix{SubgroupIndex}(undef, num_rows_conj_matrix, num_cols_conj_matrix)
        right_conj_matx =
            Matrix{SubgroupIndex}(undef, num_rows_conj_matrix, num_cols_conj_matrix)

        for i in eachindex(subgroups)
            for j in eachindex(generators)
                H = subgroups[i]
                g = generators[j]

                left_conjugated = H^g
                right_conjugated = H^(g^-1)

                for k in eachindex(subgroups)
                    if subgroups[k] == left_conjugated
                        left_conj_matx[i, j] = k
                    end

                    if subgroups[k] == right_conjugated
                        right_conj_matx[i, j] = k
                    end
                end
            end
        end

        ########## double coset stuff below here
        doubleCosetRepresentatives =
            Dict{Tuple{SubgroupIndex,SubgroupIndex,SubgroupIndex},Vector{GeneratorWord}}()

        for k in eachindex(subgroups)
            covered_subgroups =
                SubgroupIndex[h for (h, upper) in covers if upper == k]
            isempty(covered_subgroups) && continue

            K = subgroups[k]
            for h in covered_subgroups, j in covered_subgroups
                H = subgroups[h]
                J = subgroups[j]

                doubleCosetRepresentatives[(h, k, j)] = [
                    generator_word(epi_from_free_group, entry[1])
                    for entry in GAP.Globals.DoubleCosetRepsAndSizes(K, H, J)
                ]
            end
        end

        return new(
            G,
            subgroups,
            covers,
            epi_from_free_group,
            generators,
            left_conj_matx,
            right_conj_matx,
            doubleCosetRepresentatives,
        )
    end
end

"""
    double_coset_representative_words(ctx, h, k, j)
    double_coset_representative_words(ctx, H, K, J)

Return representatives for `H \\ K / J`, where `H < K` and `J < K` are cover
relations in `ctx`.

Each representative is a word in `ctx.generators`, stored as a vector of
`(generator_index, exponent)` pairs.
"""
function double_coset_representative_words(
    ctx::MackeyContext,
    h::SubgroupIndex,
    k::SubgroupIndex,
    j::SubgroupIndex,
)
    h in eachindex(ctx.subgroups) || throw(ArgumentError("h is not a subgroup index"))
    k in eachindex(ctx.subgroups) || throw(ArgumentError("k is not a subgroup index"))
    j in eachindex(ctx.subgroups) || throw(ArgumentError("j is not a subgroup index"))

    key = (h, k, j)
    haskey(ctx.doubleCosetRepresentatives, key) ||
        throw(ArgumentError(
            "double-coset representatives are stored only when (h, k) and (j, k) are cover relations",
        ))

    return ctx.doubleCosetRepresentatives[key]
end

function double_coset_representative_words(
    ctx::MackeyContext,
    H::Group,
    K::Group,
    J::Group,
)
    return double_coset_representative_words(
        ctx,
        subgroup_index(ctx, H, "H"),
        subgroup_index(ctx, K, "K"),
        subgroup_index(ctx, J, "J"),
    )
end

function generator_word(free_group_map, element::GroupElement)
    word = GAP.Globals.PreImagesRepresentative(free_group_map, element)
    external_representation = Vector{Int}(GAP.Globals.ExtRepOfObj(word))

    iseven(length(external_representation)) ||
        throw(ArgumentError("GAP returned an invalid free-group word representation"))

    result = GeneratorWord()
    for i in 1:2:length(external_representation)
        push!(
            result,
            (
                GeneratorIndex(external_representation[i]),
                Int(external_representation[i + 1]),
            ),
        )
    end

    return result
end

function subgroup_index(ctx::MackeyContext, H::Group, name::AbstractString)
    matches = findall(K -> K == H, ctx.subgroups)

    length(matches) == 1 ||
        throw(ArgumentError("$name is not uniquely represented in the subgroup list"))

    return SubgroupIndex(only(matches))
end
