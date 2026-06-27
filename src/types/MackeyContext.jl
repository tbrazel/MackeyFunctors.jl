const SubgroupIndex::DataType = Int

"""
    DoubleCosetInfo

Data attached to one representative `x` of a double coset `L \\ A / R`.

The same representative is used in two different Mackey-functor computations:
the double-coset formula uses `L^x ∩ R`, while shifts use the stabilizer
`L ∩ xRx^-1` of the point `(L, xR)` in `(G/L) x (G/R)`. GAP writes `K^g` for
`g^-1*K*g`, so the second subgroup is `L ∩ R^(x^-1)`.
"""
struct DoubleCosetInfo
    # A representative x in the ambient subgroup A for a double coset L*x*R.
    representative::GroupElement

    # Index of L^x ∩ R = x^-1*L*x ∩ R. This is the intersection that appears
    # on the restriction side of the Mackey double-coset formula.
    left_conjugate_intersection_index::SubgroupIndex

    # Index of L ∩ x*R*x^-1. For shifts this is the stabilizer of the product
    # orbit representative (L, xR); it is also the conjugate by x of the
    # subgroup stored in left_conjugate_intersection_index.
    left_intersection_conjugated_right_index::SubgroupIndex
end

"""
    MackeyContext(G)

Precompute subgroup-lattice and double-coset data for a finite GAP group `G`. This front-loads the computational effort for the group and subgroup-lattice data used when verifying the axioms for any ``G``-Mackey functor.

This method takes a group ``G`` and precomputes its list of subgroups, the list of *covers* (meaning proper subgroups ``H\\le K`` where there are no intermediate subgroups), a list of generators for the group, the relators for those generators, and matrices for how these generators act via conjugation on the list of all subgroups.

Importantly, the `MackeyContext` type also stores all the data needed to verify the double coset formula for an input Mackey functor. An important lemma is that the double coset formula can be checked along subgroups ``J\\le H \\ge K`` where each inclusion is a cover. Given a triple of subgroups as above, we have
```math
H = \\coprod_x JxK
```
and we store this as a vector of `DoubleCosetInfo` values.
"""
struct MackeyContext
    group::Group
    subgroups::Vector{Group}
    covers::Vector{Tuple{SubgroupIndex,SubgroupIndex}}
    paths::Dict{Tuple{SubgroupIndex,SubgroupIndex},Vector{Int}}
    generators::Vector{GroupElement}
    fp_isomorphism::GapObj
    generator_relations::Vector{GeneratorWord}
    generator_left_conjugation_matrix::Matrix{SubgroupIndex}
    generator_right_conjugation_matrix::Matrix{SubgroupIndex}
    double_coset_info_cache::Dict{
        Tuple{SubgroupIndex,SubgroupIndex,SubgroupIndex},
        Vector{DoubleCosetInfo},
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
        generators = Vector{GroupElement}(GAP.Globals.MinimalGeneratingSet(G))

        # Store a presentation of G on the chosen generators.  GAP's
        # Factorization uses the group object's GeneratorsOfGroup attribute,
        # which may differ from the generators stored in this context.  This
        # isomorphism lets us convert group elements into words in exactly the
        # generator coordinates used by all context matrices.
        fp_isomorphism = GAP.Globals.IsomorphismFpGroupByGenerators(
            G,
            GapObj(generators; recursive=true),
        )
        relation_words = generator_relations_from_isomorphism(fp_isomorphism)

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
        double_coset_info_cache =
            Dict{
                Tuple{SubgroupIndex,SubgroupIndex,SubgroupIndex},
                Vector{DoubleCosetInfo},
            }()

        for (h, H) in enumerate(subgroups)
            covered_subgroups =
                SubgroupIndex[j for (j, upper) in covers if upper == h]
            isempty(covered_subgroups) && continue

            for j in covered_subgroups, k in covered_subgroups
                double_coset_info_cache[(j, h, k)] =
                    _compute_double_coset_infos(subgroups, j, h, k)
            end
        end

        return new(
            G,
            subgroups,
            covers,
            paths,
            generators,
            fp_isomorphism,
            relation_words,
            left_conj_matx,
            right_conj_matx,
            double_coset_info_cache,
        )
    end
end

function Base.:(==)(a::MackeyContext, b::MackeyContext)
    return a.group == b.group &&
        a.subgroups == b.subgroups &&
        a.covers == b.covers &&
        a.paths == b.paths &&
        a.generators == b.generators &&
        a.generator_relations == b.generator_relations &&
        a.generator_left_conjugation_matrix == b.generator_left_conjugation_matrix &&
        a.generator_right_conjugation_matrix == b.generator_right_conjugation_matrix
end

function generator_word(ctx::MackeyContext, element::GroupElement)::GeneratorWord
    return generator_word_from_isomorphism(ctx.fp_isomorphism, element)
end

function generator_relations(ctx::MackeyContext)::Vector{GeneratorWord}
    return ctx.generator_relations
end

"""
    double_coset_infos(ctx, j, h, k)
    double_coset_infos(ctx, J, H, K)

Return cached data for `J \\ H / K`.

Each entry is a `DoubleCosetInfo`. The key `(j, h, k)` means the left,
ambient, and right subgroups are `ctx.subgroups[j]`, `ctx.subgroups[h]`, and
`ctx.subgroups[k]`, respectively.
"""
function double_coset_infos(
    ctx::MackeyContext,
    j::SubgroupIndex,
    h::SubgroupIndex,
    k::SubgroupIndex,
)
    checkbounds(ctx.subgroups, j)
    checkbounds(ctx.subgroups, h)
    checkbounds(ctx.subgroups, k)

    is_subgroup(ctx, j, h) ||
        throw(ArgumentError("The left subgroup must be contained in the ambient subgroup."))
    is_subgroup(ctx, k, h) ||
        throw(ArgumentError("The right subgroup must be contained in the ambient subgroup."))

    key = (j, h, k)
    if !haskey(ctx.double_coset_info_cache, key)
        ctx.double_coset_info_cache[key] =
            _compute_double_coset_infos(ctx.subgroups, j, h, k)
    end

    return ctx.double_coset_info_cache[key]
end

function double_coset_infos(
    ctx::MackeyContext,
    J::Group,
    H::Group,
    K::Group,
)
    return double_coset_infos(
        ctx,
        subgroup_index(ctx, J),
        subgroup_index(ctx, H),
        subgroup_index(ctx, K),
    )
end

double_coset_representative_data(ctx::MackeyContext, args...) =
    double_coset_infos(ctx, args...)

"""
    double_coset_representative_words(ctx, j, h, k)
    double_coset_representative_words(ctx, J, H, K)

Return only the generator words from `double_coset_representative_data`.
"""
function double_coset_representative_words(ctx::MackeyContext, args...)
    return [
        generator_word(ctx, info.representative)
        for info in double_coset_infos(ctx, args...)
    ]
end

"""
    subgroup_inclusion_index(ctx, i, j)
    subgroup_inclusion_index(ctx, (i, j))

Return the group-theoretic index ``[H[j] : H[i]]`` as an `Int`, where
`H[i]` and `H[j]` are entries of `ctx.subgroups`.
"""
function subgroup_inclusion_index(
    ctx::MackeyContext,
    i::SubgroupIndex,
    j::SubgroupIndex,
)::Int
    checkbounds(ctx.subgroups, i)
    checkbounds(ctx.subgroups, j)

    is_subgroup(ctx, i, j) ||
        throw(ArgumentError("Subgroup $i must be contained in subgroup $j."))

    return Int(GAP.Globals.Index(ctx.subgroups[j], ctx.subgroups[i]))
end

function subgroup_inclusion_index(
    ctx::MackeyContext,
    cover::Tuple{SubgroupIndex,SubgroupIndex},
)::Int
    return subgroup_inclusion_index(ctx, cover[1], cover[2])
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
                result = context.generator_left_conjugation_matrix[g, result]
            end
        else
            for j in 1:abs(n)
                result = context.generator_right_conjugation_matrix[g, result]
            end
        end
    end
    return result
end

function _compute_double_coset_infos(
    subgroups::Vector{Group},
    left_index::SubgroupIndex,
    ambient_index::SubgroupIndex,
    right_index::SubgroupIndex,
)::Vector{DoubleCosetInfo}
    left = subgroups[left_index]
    ambient = subgroups[ambient_index]
    right = subgroups[right_index]

    infos = DoubleCosetInfo[]
    for entry in GAP.Globals.DoubleCosetRepsAndSizes(ambient, left, right)
        x = entry[1]

        # Mackey formulas use L^x ∩ R. This subgroup lies inside R.
        left_conjugate_intersection =
            GAP.Globals.Intersection(left^x, right)

        # Shift decompositions use the stabilizer of (L, xR), namely
        # L ∩ xRx^-1. GAP writes xRx^-1 as R^(x^-1).
        left_intersection_conjugated_right =
            GAP.Globals.Intersection(left, right^(x^-1))

        push!(
            infos,
            DoubleCosetInfo(
                x,
                subgroup_index(subgroups, left_conjugate_intersection),
                subgroup_index(subgroups, left_intersection_conjugated_right),
            ),
        )
    end

    return infos
end

function whole_group_index(ctx::MackeyContext)::SubgroupIndex
    return subgroup_index(ctx, ctx.group)
end

# Build the double-coset data for H\G/K used to decompose
# (G/H) x (G/K) into transitive G-orbits.
function _shift_decomposition_double_coset_infos(
    ctx::MackeyContext,
    H_index::SubgroupIndex,
    K_index::SubgroupIndex,
)
    return double_coset_infos(ctx, H_index, whole_group_index(ctx), K_index)
end
