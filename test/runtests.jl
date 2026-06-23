using Test
using GAP
using AbstractAlgebra
using MackeyFunctors

function zero_mackey_functor_data(G)
    subs = MackeyFunctors.lattice_subgroups(G)
    covers = MackeyFunctors.cover_relations(subs)
    zero_map = zeros(Int64, 0, 0)

    values = IdDict{GapObj, Vector{Int64}}()
    restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}()
    transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}()
    conjugations = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}()

    for H in subs
        values[H] = Int64[]
    end

    for cover in covers
        restrictions[cover] = zero_map
        transfers[cover] = zero_map
    end

    for H in subs, g in GAP.Globals.GeneratorsOfGroup(G)
        Bool(GAP.Globals.IN(g, H)) && continue
        conjugations[(H, g)] = zero_map
    end

    return values, restrictions, transfers, conjugations
end

function c2_burnside_data()
    C2 = GAP.Globals.CyclicGroup(2)
    subs = Vector{GapObj}(GAP.Globals.AllSubgroups(C2))
    e = subs[1]
    g = GAP.Globals.GeneratorsOfGroup(C2)[1]

    values = IdDict{GapObj, Vector{Int64}}(
        e => [0],
        C2 => [0, 0]
    )
    restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [1 2]
    )
    transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [0; 1;;]
    )
    conjugations = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, g) => [1;;]
    )

    return C2, e, g, values, restrictions, transfers, conjugations
end


@testset "Abelian group maps" begin
    source_group = [0, 2, 4]   # Z ⊕ C₂ ⊕ C₄
    target_group = [3, 0, 0]   # C₃ ⊕ Z ⊕ Z
    @testset "valid map" begin
        M = [
            1  3  6
            2  0  0
            5  0  0
        ]
        @test MackeyFunctors.defines_fg_abelian_map(source_group, target_group, M)
    end

    # todo - write more tests
end

@testset "Make C2 Burnside Mackey Functor" begin
    C2, identity_subgroup, _, values, restrictions, transfers, conjugations =
        c2_burnside_data()

    M = MackeyFunctor(C2, values, restrictions, transfers, conjugations)
    e = MackeyFunctors.canonical_subgroup(M.subgroups, identity_subgroup)
    c2 = MackeyFunctors.canonical_subgroup(M.subgroups, C2)
    res_key = MackeyFunctors.gap_pair_key(M.restrictions, e, c2)
    tr_key = MackeyFunctors.gap_pair_key(M.transfers, e, c2)

    @test M.G == C2
    @test M.coefficient_ring === ZZ
    @test M.lattice isa GapObj
    @test length(M.subgroups) == 2
    @test length(M.covers) == 1
    @test all(H -> any(K -> H === K, M.subgroups), keys(M.values))
    @test M.values[e] == [0]
    @test M.values[c2] == [0, 0]
    @test res_key !== nothing
    @test tr_key !== nothing
    @test M.restrictions[res_key] == [1 2]
    @test M.transfers[tr_key] == [0; 1;;]

    F5 = GF(5)
    Zx, _ = polynomial_ring(ZZ, :x)
    M_F5 = MackeyFunctor(C2, values, restrictions, transfers, conjugations; coefficient_ring=F5)
    M_Zx = MackeyFunctor(Zx, C2, values, restrictions, transfers, conjugations)

    @test M_F5.coefficient_ring === F5
    @test M_Zx.coefficient_ring === Zx
end


@testset "Make C3 underlying free Mackey Functor" begin
    C3 = GAP.Globals.CyclicGroup(3)
    subs = Vector{GapObj}(GAP.Globals.AllSubgroups(C3))
    identity_subgroup = subs[1]
    g = GAP.Globals.GeneratorsOfGroup(C3)[1]

    c3_value = [0]
    e_value = [0,0,0]
    res = [1;1;1;;]
    tr = [1 1 1]
    conj = [0 0 1; 1 0 0; 0 1 0]

    values = IdDict{GapObj, Vector{Int64}}(
        identity_subgroup => e_value,
        C3 => c3_value
    )

    restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (identity_subgroup,C3) => res
    )

    transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (identity_subgroup,C3) => tr
    )

    conjugations = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (identity_subgroup,g) => conj
    )

    M = MackeyFunctor(C3, values, restrictions, transfers, conjugations)
    e = MackeyFunctors.canonical_subgroup(M.subgroups, identity_subgroup)
    conj_squared_key = MackeyFunctors.gap_pair_key(M.conjugations, e, g*g)
 
    @test conj_squared_key !== nothing
    @test M.conjugations[conj_squared_key] == conj * conj
end

@testset "Lattice cover data generates all comparable maps" begin
    C4 = GAP.Globals.CyclicGroup(4)
    values, restrictions, transfers, conjugations = zero_mackey_functor_data(C4)

    M = MackeyFunctor(C4, values, restrictions, transfers, conjugations)

    noncover = first(
        (H, K) for H in M.subgroups for K in M.subgroups
        if MackeyFunctors.is_proper_subgroup(H, K) &&
           !MackeyFunctors.has_cover(M.covers, H, K)
    )

    @test MackeyFunctors.gap_pair_key(M.restrictions, noncover...) !== nothing
    @test MackeyFunctors.gap_pair_key(M.transfers, noncover...) !== nothing
end

@testset "Lattice MackeyFunctor rejects invalid data" begin
    C2, e, _, values, restrictions, transfers, conjugations = c2_burnside_data()

    missing_values = copy(values)
    delete!(missing_values, C2)
    @test_throws ArgumentError MackeyFunctor(C2, missing_values, restrictions, transfers, conjugations)

    bad_restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [1; 2;;]
    )
    @test_throws DimensionMismatch MackeyFunctor(C2, values, bad_restrictions, transfers, conjugations)

    bad_transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [1; 0;;]
    )
    @test_throws ArgumentError MackeyFunctor(C2, values, restrictions, bad_transfers, conjugations)
    @test_throws ArgumentError MackeyFunctor(:not_a_ring, C2, values, restrictions, transfers, conjugations)

    C4 = GAP.Globals.CyclicGroup(4)
    zero_values, zero_restrictions, zero_transfers, zero_conjugations = zero_mackey_functor_data(C4)
    lattice_subs = MackeyFunctors.lattice_subgroups(C4)
    covers = MackeyFunctors.cover_relations(lattice_subs)
    noncover = first(
        (H, K) for H in lattice_subs for K in lattice_subs
        if MackeyFunctors.is_proper_subgroup(H, K) &&
           !MackeyFunctors.has_cover(covers, H, K)
    )

    zero_restrictions[noncover] = zeros(Int64, 0, 0)
    @test_throws ArgumentError MackeyFunctor(
        C4,
        zero_values,
        zero_restrictions,
        zero_transfers,
        zero_conjugations,
    )
end

@testset "MackeyFunctor construction examples" begin
    C2, _, _, values, restrictions, transfers, conjugations = c2_burnside_data()
    C4 = GAP.Globals.CyclicGroup(4)
    zero_values, zero_restrictions, zero_transfers, zero_conjugations =
        zero_mackey_functor_data(C4)

    @test MackeyFunctor(C2, values, restrictions, transfers, conjugations).coefficient_ring === ZZ
    @test MackeyFunctor(
        C2,
        values,
        restrictions,
        transfers,
        conjugations;
        coefficient_ring=GF(7),
    ).coefficient_ring === GF(7)
    @test length(MackeyFunctor(
        C4,
        zero_values,
        zero_restrictions,
        zero_transfers,
        zero_conjugations,
    ).subgroups) == 3
end

@testset "Mackey identity failures" begin
    C2, e, _, values, restrictions, transfers, conjugations = c2_burnside_data()

    bad_mackey_transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e, C2) => [1; 0;;]
    )
    @test_throws ArgumentError MackeyFunctor(
        C2,
        values,
        restrictions,
        bad_mackey_transfers,
        conjugations,
    )

    C3 = GAP.Globals.CyclicGroup(3)
    subs = Vector{GapObj}(GAP.Globals.AllSubgroups(C3))
    e3 = subs[1]
    g = GAP.Globals.GeneratorsOfGroup(C3)[1]

    c3_values = IdDict{GapObj, Vector{Int64}}(
        e3 => [0, 0, 0],
        C3 => [0]
    )
    c3_restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e3, C3) => [1; 1; 1;;]
    )
    c3_transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e3, C3) => [1 1 1]
    )
    bad_conjugations = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
        (e3, g) => [2 0 0; 0 1 0; 0 0 1]
    )

    @test_throws ArgumentError MackeyFunctor(
        C3,
        c3_values,
        c3_restrictions,
        c3_transfers,
        bad_conjugations,
    )
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

    relations = generator_relations(G, ctx.generators)

    @test ctx.G == G
    @test length(ctx.subgroups) == length(GAP.Globals.AllSubgroups(G))
    @test length(ctx.generators) == length(GAP.Globals.GeneratorsOfGroup(G))
    @test size(ctx.generatorLeftConjugationMatrix) ==
          (length(ctx.subgroups), length(ctx.generators))
    @test size(ctx.generatorRightConjugationMatrix) ==
          (length(ctx.subgroups), length(ctx.generators))
    @test generator_relations(G) == relations
    @test !isempty(relations)
    @test all(relations) do relation
        relation isa Vector{Tuple{Int,Int}} &&
            evaluate_word(relation) == GAP.Globals.One(G) &&
            all(relation) do (generator_index, exponent)
                1 <= generator_index <= length(ctx.generators) && exponent isa Int
            end
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

@testset "Visualizer data" begin
    C2, _, _, values, restrictions, transfers, conjugations = c2_burnside_data()
    M = MackeyFunctor(C2, values, restrictions, transfers, conjugations)
    data = visualizer_data(M)
    json = visualizer_json(M)

    @test data["coefficient_ring"] == "Integers"
    @test length(data["nodes"]) == 2
    @test length(data["covers"]) == 1
    @test length(data["restrictions"]) == 1
    @test length(data["transfers"]) == 1
    @test length(data["conjugations"]) == 2
    @test haskey(first(data["nodes"]), "label_tex")
    @test haskey(first(data["nodes"]), "value_tex")
    @test haskey(first(data["nodes"]), "orbit_tex")
    @test occursin("\"nodes\"", json)
    @test occursin("\"covers\"", json)
end
