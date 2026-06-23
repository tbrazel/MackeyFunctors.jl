using Test
using GAP
using AbstractAlgebra
using MackeyFunctors

@testset "Constant Z C2-Mackey functor" begin
    C2 = GAP.Globals.CyclicGroup(2)
    context = MackeyContext(C2)

    @test length(context.covers) == 1

    Z = free_module(ZZ,1)
    values = [Z]
    restrictions = [ModuleHomomorphism(Z,Z,identity_matrix(ZZ,1))]
    transfers = [ModuleHomomorphism(Z,Z,matrix(ZZ,[[2]]))]
    conjugations = [ModuleHomomorphism(Z,Z,identity_matrix(ZZ,1))]

    MackeyFunctor(context, values, restrictions, transfers, conjugations)
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

    @test ctx.G == G
    @test length(ctx.subgroups) == length(GAP.Globals.AllSubgroups(G))
    @test length(ctx.generators) == length(GAP.Globals.GeneratorsOfGroup(G))
    @test size(ctx.generatorLeftConjugationMatrix) ==
          (length(ctx.subgroups), length(ctx.generators))
    @test size(ctx.generatorRightConjugationMatrix) ==
          (length(ctx.subgroups), length(ctx.generators))

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
