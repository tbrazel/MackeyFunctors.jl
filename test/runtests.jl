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
            #     values = [M for i in context.subgroups]

            #     id_hom = ModuleHomomorphism(M, M, identity_matrix(R, 1))
            #     restrictions = [id_hom for i in context.covers]
            #     transfers = [R(2) * id_hom for i in context.covers]

            #     id_iso = ModuleIsomorphism(M, M, identity_matrix(R, 1))
            #     conjugations = [id_iso for i in context.generators, j in context.subgroups]

            #     @test MackeyFunctor(context, values, restrictions, transfers, conjugations) isa MackeyFunctor
            cmf = constant_mackey_functor(context, M)
            @test cmf isa MackeyFunctor
            @test coefficient_ring(cmf) == R
        end
    end
end

@testset "Constant Mackey functors for symmetric groups" begin
    for n in 2:5
        S_n = GAP.Globals.SymmetricGroup(n)
        context = MackeyContext(S_n)

        @test length(context.subgroups) == length(GAP.Globals.AllSubgroups(S_n))

        for R in [ZZ, QQ, GF(2), GF(67)]
            M = free_module(R, 1)
            @test constant_mackey_functor(context, M) isa MackeyFunctor
        end
    end
end

@testset "Burnside Mackey functors for Cp" begin
    for p in [2, 3, 5, 7, 11, 67]
        Cp = GAP.Globals.CyclicGroup(p)
        context = MackeyContext(Cp)

        @test length(context.subgroups) == 2
        @test length(context.covers) == 1

        for R in [ZZ, QQ, GF(p), GF(1000000007)]
            Ae = free_module(R, 1)
            ACp = free_module(R, 2)
            val = [Ae, ACp]

            res = [hom(ACp, Ae, R[1; p])]
            tr = [hom(Ae, ACp, R[0 1;])]

            conj_e = ModuleIsomorphism(Ae, Ae, R[1;])
            conj_Cp = ModuleIsomorphism(ACp, ACp, identity_matrix(R, 2))
            conj = Generic.ModuleIsomorphism[conj_e conj_Cp;]

            @test MackeyFunctor(context, val, res, tr, conj) isa MackeyFunctor
        end
    end
end

# eventually we will have a constructor for these,
# but putting this here for us to compare against
@testset "Free Mackey functors for C4" begin
    C4 = GAP.Globals.CyclicGroup(4)
    context = MackeyContext(C4)

    for R in [ZZ, QQ, GF(2), GF(67)]
        # A4 = Burnside Mackey functor
        A4_1 = free_module(R, 1)
        A4_2 = free_module(R, 2)
        A4_4 = free_module(R, 3)
        A4_val = [A4_1, A4_2, A4_4]

        A4_res = [hom(A4_2, A4_1, R[1; 2]),
            hom(A4_4, A4_2, R[1 0; 2 0; 0 2])]
        A4_tr = [hom(A4_1, A4_2, R[0 1;]),
            hom(A4_2, A4_4, R[0 1 0; 0 0 1])]
        # all conjugations trivial in A2
        A4_conj = [MackeyFunctors.identity_isomorphism(A4_val[j])
            for j in eachindex(context.generators), j in eachindex(A4_val)]

        A4 = MackeyFunctor(context, A4_val, A4_res, A4_tr, A4_conj)
        @test A4 isa MackeyFunctor

        # A2 = free at level C2
        A2_1 = free_module(R, 2)
        A2_2 = free_module(R, 4)
        A2_4 = free_module(R, 2)
        A2_val = [A2_1, A2_2, A2_4]

        A2_res = [hom(A2_2, A2_1, R[1 0; 2 0; 0 1; 0 2]),
            hom(A2_4, A2_2, R[1 0 1 0; 0 1 0 1])]
        A2_tr = [hom(A2_1, A2_2, R[0 1 0 0; 0 0 0 1]),
            hom(A2_2, A2_4, R[1 0; 0 1; 1 0; 0 1])]
        A2_conj = Generic.ModuleIsomorphism[ #=
            =# ModuleIsomorphism(A2_1, A2_1,  R[0 1; 1 0]) #=
            =# ModuleIsomorphism(A2_2, A2_2, R[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]) #=
            =# MackeyFunctors.identity_isomorphism(A2_4); #=
            =# MackeyFunctors.identity_isomorphism(A2_1) #=
            =# MackeyFunctors.identity_isomorphism(A2_2) #=
            =# MackeyFunctors.identity_isomorphism(A2_4)]

        A2 = MackeyFunctor(context, A2_val, A2_res, A2_tr, A2_conj)
        @test A2 isa MackeyFunctor

        # A1 = free at underlying level
        A1_1 = free_module(R, 4)
        A1_2 = free_module(R, 2)
        A1_4 = free_module(R, 1)
        A1_val = [A1_1, A1_2, A1_4]

        A1_res = [hom(A1_2, A1_1, R[1 0 1 0; 0 1 0 1]),
            hom(A1_4, A1_2, R[1 1;])]
        A1_tr = [hom(A1_1, A1_2, R[1 0; 0 1; 1 0; 0 1]),
            hom(A1_2, A1_4, R[1; 1])]
        A1_conj = Generic.ModuleIsomorphism[ #=
            =# ModuleIsomorphism(A1_1, A1_1, R[0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0]) #=
            =# ModuleIsomorphism(A1_2, A1_2, R[0 1; 1 0]) #=
            =# MackeyFunctors.identity_isomorphism(A1_4); #=
            =# ModuleIsomorphism(A1_1, A1_1,  R[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]) #=
            =# MackeyFunctors.identity_isomorphism(A1_2) #=
            =# MackeyFunctors.identity_isomorphism(A1_4)]

        A1 = MackeyFunctor(context, A1_val, A1_res, A1_tr, A1_conj)
        @test A1 isa MackeyFunctor

        # TODO: uncomment when free constructor is done
        # check A1 is free_mackey_functor(context, 1)
        # check A2 is free_mackey_functor(context, 2)
        # check A4 is free_mackey_functor(context, 3)
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
    conj = MackeyFunctors.conjugation(s3_mackey_functor, conjugation_word, 1)
    @test all(gens(M)) do x
        conj(x) == A_iso(B_iso(x))
    end
end

@testset "Shifts of Mackey functors" begin
    C4 = GAP.Globals.CyclicGroup(4)
    c4_context = MackeyContext(C4)
    c4_constant = constant_mackey_functor(c4_context, ZZ)
    c4_trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, c4_context.subgroups)
    c4_whole = findfirst(
        H -> Int(GAP.Globals.Size(H)) == Int(GAP.Globals.Size(C4)),
        c4_context.subgroups,
    )

    shifted_by_trivial = shift(c4_constant, c4_trivial)
    @test shifted_by_trivial isa MackeyFunctor
    @test shift(c4_constant, c4_trivial; verify=false) isa MackeyFunctor

    for (K_index, K) in enumerate(c4_context.subgroups)
        expected_rank = Int(GAP.Globals.Index(C4, K))
        @test rank(MackeyFunctors.value(shifted_by_trivial, K_index)) == expected_rank
    end

    shifted_by_whole = shift(c4_constant, c4_whole)
    @test shifted_by_whole isa MackeyFunctor
    @test all(eachindex(c4_context.subgroups)) do K_index
        rank(MackeyFunctors.value(shifted_by_whole, K_index)) == 1
    end
    @test all(eachindex(c4_context.covers)) do cover_index
        matrix(shifted_by_whole.cover_restrictions[cover_index]) ==
            matrix(c4_constant.cover_restrictions[cover_index]) &&
            matrix(shifted_by_whole.cover_transfers[cover_index]) ==
            matrix(c4_constant.cover_transfers[cover_index])
    end

    S3 = GAP.Globals.SymmetricGroup(3)
    s3_context = MackeyContext(S3)
    s3_constant = constant_mackey_functor(s3_context, ZZ)
    subgroup_of_order_two =
        findfirst(H -> Int(GAP.Globals.Size(H)) == 2, s3_context.subgroups)

    shifted_s3 = shift(s3_constant, subgroup_of_order_two)
    @test shifted_s3 isa MackeyFunctor

    for (K_index, K) in enumerate(s3_context.subgroups)
        expected_rank = length(
            GAP.Globals.DoubleCosetRepsAndSizes(
                S3,
                s3_context.subgroups[subgroup_of_order_two],
                K,
            ),
        )
        @test rank(MackeyFunctors.value(shifted_s3, K_index)) == expected_rank
    end

    s3_burnside = burnside_mackey_functor(s3_context, ZZ)
    @test shift(s3_burnside, subgroup_of_order_two) isa MackeyFunctor
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

    @test Set(keys(ctx.double_coset_formulae)) == expected_double_coset_keys

    for ((j, h, k), representatives) in ctx.double_coset_formulae
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

    @test !haskey(ctx.double_coset_formulae, (trivial, whole, trivial))
    @test_throws ArgumentError double_coset_representative_data(ctx, trivial, whole, trivial)
    @test_throws ArgumentError double_coset_representative_words(ctx, trivial, whole, trivial)

    j, h = first(ctx.covers)
    @test haskey(ctx.double_coset_formulae, (j, h, j))

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
