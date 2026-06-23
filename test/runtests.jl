using Test
using GAP
using AbstractAlgebra
using MackeyFunctors

@testset "Constant Mackey functors for cyclic 2-groups" begin
    for k in 1:5
        C2k = GAP.Globals.CyclicGroup(2^k)
        context = MackeyContext(C2k)

        @test length(context.subgroups) == k+1
        @test length(context.covers) == k

        for R in [ZZ, QQ, GF(2), GF(67)]
            M = free_module(R, 1)
            values = [M for i in context.subgroups]

            id_hom = ModuleHomomorphism(M, M, identity_matrix(R, 1))
            restrictions = [id_hom for i in context.covers]
            transfers = [R(2) * id_hom for i in context.covers]

            id_iso = ModuleIsomorphism(M, M, identity_matrix(R, 1))
            conjugations = [id_iso for i in context.generators, j in context.subgroups]

            MackeyFunctor(context, values, restrictions, transfers, conjugations)
        end
    end
end

@testset "Module map composition order" begin
    M = free_module(ZZ, 2)
    A_hom = ModuleHomomorphism(M, M, matrix(ZZ, [1 1; 0 1]))
    B_hom = ModuleHomomorphism(M, M, matrix(ZZ, [1 0; 1 1]))
    A_iso = ModuleIsomorphism(M, M, matrix(ZZ, [1 1; 0 1]))
    B_iso = ModuleIsomorphism(M, M, matrix(ZZ, [1 0; 1 1]))
    id_hom = ModuleHomomorphism(M, M, identity_matrix(ZZ, 2))
    id_iso = ModuleIsomorphism(M, M, identity_matrix(ZZ, 2))

    C4 = GAP.Globals.CyclicGroup(4)
    c4_context = MackeyContext(C4)
    trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, c4_context.subgroups)
    whole = findfirst(
        H -> Int(GAP.Globals.Size(H)) == Int(GAP.Globals.Size(C4)),
        c4_context.subgroups,
    )
    path = c4_context.paths[(trivial, whole)]
    @test length(path) == 2

    values = [M for _ in c4_context.subgroups]
    restrictions = [
        i == path[1] ? A_hom : i == path[2] ? B_hom : id_hom
        for i in eachindex(c4_context.covers)
    ]
    transfers = copy(restrictions)
    conjugations = [id_iso for _ in c4_context.generators, _ in c4_context.subgroups]
    c4_mackey_functor = MackeyFunctor(
        c4_context,
        values,
        restrictions,
        transfers,
        conjugations,
        false,
    )

    res = MackeyFunctors.restriction(c4_mackey_functor, trivial, whole)
    tr = MackeyFunctors.transfer(c4_mackey_functor, trivial, whole)

    @test all(gens(M)) do x
        res(x) == restrictions[path[1]](restrictions[path[2]](x))
    end
    @test all(gens(M)) do x
        tr(x) == transfers[path[2]](transfers[path[1]](x))
    end

    S3 = GAP.Globals.SymmetricGroup(3)
    s3_context = MackeyContext(S3)
    s3_restrictions = [id_hom for _ in s3_context.covers]
    s3_transfers = [id_hom for _ in s3_context.covers]
    s3_conjugations = [
        i == 1 ? A_iso : i == 2 ? B_iso : id_iso
        for i in eachindex(s3_context.generators), _ in s3_context.subgroups
    ]
    s3_mackey_functor = MackeyFunctor(
        s3_context,
        [M for _ in s3_context.subgroups],
        s3_restrictions,
        s3_transfers,
        s3_conjugations,
        false,
    )

    conjugation_word = [(1, 1), (2, 1)]
    conj = MackeyFunctors.conjugation(s3_mackey_functor, 1, conjugation_word)
    @test all(gens(M)) do x
        conj(x) == A_iso(B_iso(x))
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
    @test length(ctx.subgroups) == length(GAP.Globals.AllSubgroups(G))
    @test length(ctx.generators) == length(GAP.Globals.GeneratorsOfGroup(G))

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

    @test Set(keys(ctx.doubleCosetRepresentatives)) == expected_double_coset_keys

    for ((j, h, k), representatives) in ctx.doubleCosetRepresentatives
        @test (j, h) in cover_pairs
        @test (k, h) in cover_pairs
        @test all(representatives) do (word, intersection_index)
            x = evaluate_word(word)
            expected_intersection =
                GAP.Globals.Intersection(ctx.subgroups[j]^x, ctx.subgroups[k])

            1 <= intersection_index <= length(ctx.subgroups) &&
                ctx.subgroups[intersection_index] == expected_intersection &&
                Bool(GAP.Globals.IN(x, ctx.subgroups[h])) &&
                all(word) do (generator_index, exponent)
                    1 <= generator_index <= length(ctx.generators) && exponent isa Int
                end
        end
        @test all(first.(representatives)) do word
            all(word) do (generator_index, exponent)
                1 <= generator_index <= length(ctx.generators) && exponent isa Int
            end
        end
    end

    trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, ctx.subgroups)
    whole = findfirst(H -> Int(GAP.Globals.Size(H)) == Int(GAP.Globals.Size(G)), ctx.subgroups)

    @test !haskey(ctx.doubleCosetRepresentatives, (trivial, whole, trivial))
    @test_throws ArgumentError double_coset_representative_data(ctx, trivial, whole, trivial)
    @test_throws ArgumentError double_coset_representative_words(ctx, trivial, whole, trivial)

    j, h = first(ctx.covers)
    @test haskey(ctx.doubleCosetRepresentatives, (j, h, j))

    representative_data = double_coset_representative_data(ctx, j, h, j)
    @test double_coset_representative_data(
        ctx,
        ctx.subgroups[j],
        ctx.subgroups[h],
        ctx.subgroups[j],
    ) == representative_data

    words = double_coset_representative_words(ctx, j, h, j)
    @test words == first.(representative_data)

    gap_representatives = [
        entry[1]
        for entry in GAP.Globals.DoubleCosetRepsAndSizes(
            ctx.subgroups[h],
            ctx.subgroups[j],
            ctx.subgroups[j],
        )
    ]
    @test length(words) == length(gap_representatives)
    @test all(zip(representative_data, gap_representatives)) do ((word, intersection_index), representative)
        evaluate_word(word) == representative &&
            ctx.subgroups[intersection_index] ==
            GAP.Globals.Intersection(ctx.subgroups[j]^representative, ctx.subgroups[j])
    end
end

# @testset "Visualizer data" begin
#     C2, _, _, values, restrictions, transfers, conjugations = c2_burnside_data()
#     M = MackeyFunctor(C2, values, restrictions, transfers, conjugations)
#     data = visualizer_data(M)
#     json = visualizer_json(M)

#     @test data["coefficient_ring"] == "Integers"
#     @test length(data["nodes"]) == 2
#     @test length(data["covers"]) == 1
#     @test length(data["restrictions"]) == 1
#     @test length(data["transfers"]) == 1
#     @test length(data["conjugations"]) == 2
#     @test haskey(first(data["nodes"]), "label_tex")
#     @test haskey(first(data["nodes"]), "value_tex")
#     @test haskey(first(data["nodes"]), "orbit_tex")
#     @test occursin("\"nodes\"", json)
#     @test occursin("\"covers\"", json)
# end
