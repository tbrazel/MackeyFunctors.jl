@testset "MackeyContext for the trivial group" begin
    G = GAP.Globals.CyclicGroup(1)
    ctx = MackeyContext(G)

    @test ctx.group == G
    @test length(ctx.subgroups) == 1
    @test isempty(ctx.covers)
    @test ctx.paths == Dict((1, 1) => Int[])
    @test isempty(ctx.double_coset_info_cache)
    @test all(ctx.generator_relations) do word
        isempty(word)
    end
end

@testset "MackeyContext" begin
    G = GAP.Globals.SymmetricGroup(3)
    ctx = MackeyContext(G)

    function evaluate_word(word)
        result = GAP.Globals.One(G)
        for (generator_index, exponent) in word
            result *= ctx.generators[generator_index]^exponent
        end
        return result
    end

    @test ctx.group == G
    @test ctx == MackeyContext(G)
    @test ctx != MackeyContext(GAP.Globals.CyclicGroup(2))
    @test length(ctx.subgroups) == length(GAP.Globals.AllSubgroups(G))
    @test length(ctx.generators) == length(GAP.Globals.MinimalGeneratingSet(G))
    @test all(ctx.generators) do g
        all(entry -> 1 <= entry[1] <= length(ctx.generators), MackeyFunctors.generator_word(ctx, g))
    end
    @test MackeyFunctors.generator_relations(ctx) === ctx.generator_relations
    @test all(ctx.generator_relations) do word
        evaluate_word(word) == GAP.Globals.One(G)
    end

    for (i, j) in ctx.covers
        H = ctx.subgroups[i]
        K = ctx.subgroups[j]

        @test i != j
        @test Bool(GAP.Globals.IsSubgroup(K, H))
        @test !any(eachindex(ctx.subgroups)) do l
            L = ctx.subgroups[l]
            l != i &&
                l != j &&
                Bool(GAP.Globals.IsSubgroup(L, H)) &&
                Bool(GAP.Globals.IsSubgroup(K, L))
        end
    end

    cover_pairs = Set(ctx.covers)
    expected_double_coset_keys = Set{Tuple{Int,Int,Int}}()
    for (j, h) in ctx.covers, (k, upper) in ctx.covers
        h == upper || continue
        push!(expected_double_coset_keys, (j, h, k))
    end

    @test Set(keys(ctx.double_coset_info_cache)) == expected_double_coset_keys

    for ((j, h, k), representatives) in ctx.double_coset_info_cache
        @test (j, h) in cover_pairs
        @test (k, h) in cover_pairs
        @test all(representatives) do info
            x = info.representative
            expected_left_conjugate_intersection =
                GAP.Globals.Intersection(ctx.subgroups[j]^x, ctx.subgroups[k])
            expected_left_intersection_conjugated_right =
                GAP.Globals.Intersection(ctx.subgroups[j], ctx.subgroups[k]^(x^-1))

            1 <= info.left_conjugate_intersection_index <= length(ctx.subgroups) &&
                1 <= info.left_intersection_conjugated_right_index <= length(ctx.subgroups) &&
                ctx.subgroups[info.left_conjugate_intersection_index] ==
                expected_left_conjugate_intersection &&
                ctx.subgroups[info.left_intersection_conjugated_right_index] ==
                expected_left_intersection_conjugated_right &&
                Bool(GAP.Globals.IN(x, ctx.subgroups[h]))
        end
        @test all(MackeyFunctors.double_coset_representative_words(ctx, j, h, k)) do word
            all(word) do (generator_index, exponent)
                1 <= generator_index <= length(ctx.generators) && exponent isa Int
            end
        end
    end

    trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, ctx.subgroups)
    whole = findfirst(H -> Int(GAP.Globals.Size(H)) == Int(GAP.Globals.Size(G)), ctx.subgroups)

    @test !haskey(ctx.double_coset_info_cache, (trivial, whole, trivial))
    general_representative_data =
        MackeyFunctors.double_coset_representative_data(ctx, trivial, whole, trivial)
    @test haskey(ctx.double_coset_info_cache, (trivial, whole, trivial))
    @test length(general_representative_data) == length(
        GAP.Globals.DoubleCosetRepsAndSizes(
            ctx.subgroups[whole],
            ctx.subgroups[trivial],
            ctx.subgroups[trivial],
        ),
    )

    j, h = first(ctx.covers)
    @test haskey(ctx.double_coset_info_cache, (j, h, j))

    representative_data = MackeyFunctors.double_coset_representative_data(ctx, j, h, j)
    @test MackeyFunctors.double_coset_representative_data(
        ctx,
        ctx.subgroups[j],
        ctx.subgroups[h],
        ctx.subgroups[j],
    ) == representative_data

    words = MackeyFunctors.double_coset_representative_words(ctx, j, h, j)

    gap_representatives = [
        entry[1]
        for entry in GAP.Globals.DoubleCosetRepsAndSizes(
            ctx.subgroups[h],
            ctx.subgroups[j],
            ctx.subgroups[j],
        )
    ]
    @test length(words) == length(gap_representatives)
    @test all(zip(words, gap_representatives)) do (word, representative)
        evaluate_word(word) == representative
    end
    @test all(zip(representative_data, gap_representatives)) do (info, representative)
        info.representative == representative &&
            ctx.subgroups[info.left_conjugate_intersection_index] ==
            GAP.Globals.Intersection(ctx.subgroups[j]^representative, ctx.subgroups[j]) &&
            ctx.subgroups[info.left_intersection_conjugated_right_index] ==
            GAP.Globals.Intersection(ctx.subgroups[j], ctx.subgroups[j]^(representative^-1))
    end
end
