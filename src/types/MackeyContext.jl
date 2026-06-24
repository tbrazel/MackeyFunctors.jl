const SubgroupIndex::DataType = Int
const DoubleCosetFormulaTerm = Tuple{GeneratorWord,SubgroupIndex}

"""
    MackeyContext(G)

Precompute subgroup-lattice and double-coset data for a finite GAP group `G`. This front-loads the computational effort of aspects of the group and its subgroup lattice which will be used throughout verifying the axioms for any ``G``-Mackey functor.

This method inputs a group ``G``, and precomputes its list of subgroups, the list of *covers* (meaning proper subgroups ``H\\le K`` where there are no intermediate subgroups), a list of generators for the group, and matrices for how these generators act via conjugation on the list of all subgroups.

Importantly, the `MackeyContext` type also stores all the data needed to verify the double coset formula for an inputted Mackey functor. An important lemma is that the double coset formula can be checked along subgroups ``J\\le H \\ge K`` where each inclusion is a cover. Given a triple of subgroups as above, we have
```math
H = \\coprod_x JxK
```
and we store this in a dictionary as a vector of entries ``(J^x\\cap K, x)``.
"""
struct MackeyContext
    group::Group
    subgroups::Vector{Group}
    covers::Vector{Tuple{SubgroupIndex,SubgroupIndex}}
    paths::Dict{Tuple{SubgroupIndex,SubgroupIndex},Vector{Int}}
    generators::Vector{GroupElement}
    generatorLeftConjugationMatrix::Matrix{SubgroupIndex}
    generatorRightConjugationMatrix::Matrix{SubgroupIndex}
    doubleCosetRepresentatives::Dict{
        Tuple{SubgroupIndex,SubgroupIndex,SubgroupIndex},
        Vector{DoubleCosetFormulaTerm},
    }


    # Build a MackeyContext from a group G
    function MackeyContext(G::Group)::MackeyContext
        GAP.Globals.IsGroup(G) || throw(ArgumentError("Input G must be a GAP group"))

        # get subgroups
        subgroups = Vector{Group}(GAP.Globals.AllSubgroups(G))

        # Build list of covers
        covers = Tuple{SubgroupIndex,SubgroupIndex}[]

        for (i, H) in enumerate(subgroups), (j, K) in enumerate(subgroups)
            # GAP.Globals.IsSubgroup(K, H) means H <= K
            if i == j || !Bool(GAP.Globals.IsSubgroup(K, H))
                continue
            end

            has_intermediate = any(enumerate(subgroups)) do (l, L)
                H_properly_contained_in_L =
                    l != i && Bool(GAP.Globals.IsSubgroup(L, H))

                L_properly_contained_in_K =
                    l != j && Bool(GAP.Globals.IsSubgroup(K, L))

                H_properly_contained_in_L && L_properly_contained_in_K
            end

            has_intermediate || push!(covers, (i, j))
        end

        paths = Dict{Tuple{SubgroupIndex,SubgroupIndex},Vector{Int}}()

        for i in eachindex(subgroups)
            push!(paths,
                (i, i) => Int[]
            )
        end

        outgoing_covers = [Tuple{Int,SubgroupIndex}[] for _ in eachindex(subgroups)]
        for (cover_index, (i, j)) in enumerate(covers)
            push!(outgoing_covers[i], (cover_index, j))
        end

        queue = collect(keys(paths))
        queue_head = 1
        while queue_head <= length(queue)
            source, target = queue[queue_head]
            queue_head += 1

            path = paths[(source, target)]
            for (cover_index, next_target) in outgoing_covers[target]
                new_key = (source, next_target)
                new_path_length = length(path) + 1
                if !haskey(paths, new_key) || new_path_length < length(paths[new_key])
                    paths[new_key] = vcat(path, cover_index)
                    push!(queue, new_key)
                end
            end
        end

        # Get list of generators
        generators = Vector{GroupElement}(GAP.Globals.GeneratorsOfGroup(G))

        # Get epimorphism from free group
        # epi_from_free_group = GAP.Globals.EpimorphismFromFreeGroup(G)

        # Build conjugation matrices
        num_rows_conj_matrix = length(generators)
        num_cols_conj_matrix = length(subgroups)
        left_conj_matx =
            Matrix{SubgroupIndex}(undef, num_rows_conj_matrix, num_cols_conj_matrix)
        right_conj_matx =
            Matrix{SubgroupIndex}(undef, num_rows_conj_matrix, num_cols_conj_matrix)

        for i in eachindex(subgroups)
            for j in eachindex(generators)
                H = subgroups[i]
                g = generators[j]

                left_conjugated = H^(g^-1)
                right_conjugated = H^(g)

                for k in eachindex(subgroups)
                    if subgroups[k] == left_conjugated
                        left_conj_matx[j, i] = k
                    end

                    if subgroups[k] == right_conjugated
                        right_conj_matx[j, i] = k
                    end
                end
            end
        end

        ########## double coset stuff below here
        doubleCosetRepresentatives =
            Dict{
                Tuple{SubgroupIndex,SubgroupIndex,SubgroupIndex},
                Vector{DoubleCosetFormulaTerm},
            }()

        for (h, H) in enumerate(subgroups)
            covered_subgroups =
                SubgroupIndex[j for (j, upper) in covers if upper == h]
            isempty(covered_subgroups) && continue

            for j in covered_subgroups, k in covered_subgroups
                J = subgroups[j]
                K = subgroups[k]

                doubleCosetRepresentatives[(j, h, k)] = [
                    begin
                        x = entry[1]
                        intersection = GAP.Globals.Intersection(J^x, K)
                        (
                            generator_word(G, x),
                            subgroup_index(subgroups, intersection),
                        )
                    end
                    for entry in GAP.Globals.DoubleCosetRepsAndSizes(H, J, K)
                ]
            end
        end

        return new(
            G,
            subgroups,
            covers,
            paths,
            generators,
            left_conj_matx,
            right_conj_matx,
            doubleCosetRepresentatives,
        )
    end
end

"""
    double_coset_representative_data(ctx, j, h, k)
    double_coset_representative_data(ctx, J, H, K)

Return representatives for `J \\ H / K`, where `J < H` and `K < H` are cover
relations in `ctx`.

Each representative is returned as `(word, intersection_index)`, where
`intersection_index` points to `J^x` intersected with `K` in `ctx.subgroups`
for the represented group element `x`.
"""
function double_coset_representative_data(
    ctx::MackeyContext,
    j::SubgroupIndex,
    h::SubgroupIndex,
    k::SubgroupIndex,
)
    checkbounds(ctx.subgroups, j)
    checkbounds(ctx.subgroups, h)
    checkbounds(ctx.subgroups, k)

    key = (j, h, k)
    haskey(ctx.doubleCosetRepresentatives, key) ||
        throw(ArgumentError(
            "double-coset representatives are stored only when (j, h) and (k, h) are cover relations",
        ))

    return ctx.doubleCosetRepresentatives[key]
end

function double_coset_representative_data(
    ctx::MackeyContext,
    J::Group,
    H::Group,
    K::Group,
)
    return double_coset_representative_data(
        ctx,
        subgroup_index(ctx, J),
        subgroup_index(ctx, H),
        subgroup_index(ctx, K),
    )
end

"""
    double_coset_representative_words(ctx, j, h, k)
    double_coset_representative_words(ctx, J, H, K)

Return only the generator words from `double_coset_representative_data`.
"""
function double_coset_representative_words(ctx::MackeyContext, args...)
    return first.(double_coset_representative_data(ctx, args...))
end

function subgroup_index(subgroups::Vector{Group}, H::Group)::SubgroupIndex
    index = findfirst(K -> K == H, subgroups)
    index === nothing &&
        throw(ArgumentError("Subgroup is not represented in the subgroup list"))

    return index
end

function is_subgroup(ctx::MackeyContext, i::SubgroupIndex, j::SubgroupIndex)
    haskey(ctx.paths, (i, j))
end

function subgroup_index(ctx::MackeyContext, H::Group)
    subgroup_index(ctx.subgroups, H)
end

function conjugate_subgroup_by_word(context::MackeyContext, i::SubgroupIndex, w::GeneratorWord)::SubgroupIndex
    result = i
    for (g, n) in reverse(w)
        if n>0
            for j in 1:n
                result = context.generatorLeftConjugationMatrix[g, result]
            end
        else
            for j in 1:abs(n)
                result = context.generatorRightConjugationMatrix[g, result]
            end
        end
    end
    return result
end
