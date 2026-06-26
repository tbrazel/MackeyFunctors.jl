@testset "Constant Mackey functors for cyclic 2-groups" begin
    for k in 0:5
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

@testset "Mackey functor value base rings" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))

    M_ZZ = free_module(ZZ, 1)
    M_QQ = free_module(QQ, 1)
    id_ZZ_hom = ModuleHomomorphism(M_ZZ, M_ZZ, identity_matrix(ZZ, 1))
    id_ZZ_iso = ModuleIsomorphism(M_ZZ, M_ZZ, identity_matrix(ZZ, 1))

    values = AbstractAlgebra.FPModule[M_ZZ, M_QQ]
    restrictions = Generic.ModuleHomomorphism[id_ZZ_hom for _ in context.covers]
    transfers = Generic.ModuleHomomorphism[id_ZZ_hom for _ in context.covers]
    conjugations = Generic.ModuleIsomorphism[
        id_ZZ_iso for _ in context.generators, _ in context.subgroups
    ]

    @test_throws ArgumentError MackeyFunctor(
        context,
        values,
        restrictions,
        transfers,
        conjugations,
    )
    @test MackeyFunctor(
        context,
        values,
        restrictions,
        transfers,
        conjugations;
        verify=false,
    ) isa MackeyFunctor
end

@testset "Constant Mackey functors for symmetric groups" begin
    for n in 1:5
        S_n = GAP.Globals.SymmetricGroup(n)
        context = MackeyContext(S_n)

        @test length(context.subgroups) == length(GAP.Globals.AllSubgroups(S_n))

        for R in [ZZ, QQ, GF(2), GF(67)]
            M = free_module(R, 1)
            @test constant_mackey_functor(context, M) isa MackeyFunctor
        end
    end
end

@testset "Burnside Mackey functors" begin
    # Test that the constructor works
    for n in 1:4, R in [ZZ, GF(2)]
        G = GAP.Globals.SymmetricGroup(n)
        @test burnside_mackey_functor(G, R) isa MackeyFunctor
    end

    # Test that the conjugation action is nontrivial in some cases.
    # GAP search over SmallGroup(order, id) finds the first example at
    # SmallGroup(8, 3), the dihedral group D8.  It has a V4 subgroup H whose
    # normalizer is all of D8, and conjugation by an element of N_G(H) moves two
    # subgroups of H that are not conjugate inside H.
    G = GAP.Globals.SmallGroup(8, 3)
    @test string(GAP.Globals.StructureDescription(G)) == "D8"

    context = MackeyContext(G)
    H_index = findfirst(eachindex(context.subgroups)) do i
        H = context.subgroups[i]
        Int(GAP.Globals.Size(H)) == 4 &&
            string(GAP.Globals.StructureDescription(H)) == "C2 x C2" &&
            Int(GAP.Globals.Size(GAP.Globals.Normalizer(G, H))) == 8
    end
    @test H_index !== nothing

    H = context.subgroups[H_index]
    normalizer = GAP.Globals.Normalizer(G, H)
    conjugacy_classes = collect(GAP.Globals.ConjugacyClassesSubgroups(H))
    witness = nothing
    for g in GAP.Globals.Elements(normalizer), (source_index, conjugacy_class) in enumerate(conjugacy_classes)
        K = GAP.Globals.Representative(conjugacy_class)
        conjugated_K = K^(g^-1)
        target_index = findfirst(
            ==(GAP.Globals.ConjugacyClassSubgroups(H, conjugated_K)),
            conjugacy_classes,
        )

        if target_index !== nothing && target_index != source_index
            witness = (g, source_index, target_index)
            break
        end
    end
    @test witness !== nothing

    g, source_index, target_index = witness
    burnside = burnside_mackey_functor(context, ZZ)
    conjugation_action = conjugation(burnside, g, H_index)
    AH = value(burnside, H_index)

    @test conjugation_action(gens(AH)[source_index]) == gens(AH)[target_index]
    @test matrix(conjugation_action) != identity_matrix(ZZ, rank(AH))
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

@testset "Free Mackey functors for C4" begin
    C4 = GAP.Globals.CyclicGroup(4)
    context = MackeyContext(C4)

    for R in [ZZ, QQ, GF(2), GF(67)]
        # Hardcoded free Mackey functors for C4
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
            =# MackeyFunctors.identity_isomorphism(A1_4)]

        A1 = MackeyFunctor(context, A1_val, A1_res, A1_tr, A1_conj)
        @test A1 isa MackeyFunctor

        # TODO: uncomment when free constructor is done and when == is implemented??
        # check A1 is free_mackey_functor(context, 1)
        # check A2 is free_mackey_functor(context, 2)
        # check A4 is free_mackey_functor(context, 3)

        # checking ranks of Burnside shifts
        A = burnside_mackey_functor(C4, R)
        A1_free = free_mackey_functor(context, 1, R)

        A1_sh = shift(A,1)
        @test rank(A1_sh.values[1]) == 4
        @test rank(A1_sh.values[2]) == 2
        @test rank(A1_sh.values[3]) == 1
        @test all(eachindex(context.subgroups)) do i
            rank(A1_free.values[i]) == rank(A1_sh.values[i])
        end

        A2_sh = shift(A,2)
        @test rank(A2_sh.values[1]) == 2
        @test rank(A2_sh.values[2]) == 4
        @test rank(A2_sh.values[3]) == 2
    end
end

@testset "Shifts of Mackey functors for C4, S3" begin
    C1 = GAP.Globals.CyclicGroup(1)
    c1_context = MackeyContext(C1)
    c1_constant = constant_mackey_functor(c1_context, ZZ)
    shifted_c1 = shift(c1_constant, 1)

    @test shifted_c1 isa MackeyFunctor
    @test length(c1_context.covers) == 0
    @test length(shifted_c1.cover_restrictions) == 0
    @test length(shifted_c1.cover_transfers) == 0
    @test rank(MackeyFunctors.value(shifted_c1, 1)) == 1

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

    C2 = GAP.Globals.CyclicGroup(2)
    c2_context = MackeyContext(C2)
    c2_trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, c2_context.subgroups)
    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    Z2, = quo(F1, twoF1)
    c2_z2_constant = constant_mackey_functor(c2_context, Z2)
    shifted_z2 = shift(c2_z2_constant, c2_trivial; verify=false)
    shifted_z2_at_trivial = MackeyFunctors.value(shifted_z2, c2_trivial)
    @test AbstractAlgebra.invariant_factors(shifted_z2_at_trivial) == BigInt[2, 2]
    @test Set([[relation[1, i] for i in 1:ncols(relation)] for relation in relations(shifted_z2_at_trivial)]) ==
        Set([[ZZ(2), ZZ(0)], [ZZ(0), ZZ(2)]])

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

@testset "permutation_module" begin
    for n in [1, 4], R in [ZZ, GF(5)]
        G = GAP.Globals.SymmetricGroup(n)
        mc = MackeyContext(G)
        gm = permutation_module(mc, R)
        @test gm isa GModule
    end
end

@testset "fixedpoint_mackey_functor" begin
    for n in [1, 4], R in [ZZ, GF(5)]
        G = GAP.Globals.SymmetricGroup(n)
        mc = MackeyContext(G)
        gm = permutation_module(mc, R)
        mf = fixedpoint_mackey_functor(gm)
        @test mf isa MackeyFunctor
    end
end

@testset "direct sums of Mackey functors" begin
    function direct_sum_test_matrix(p, q)
        R = base_ring(domain(p))
        result = zero_matrix(
            R,
            ngens(domain(p)) + ngens(domain(q)),
            ngens(codomain(p)) + ngens(codomain(q)),
        )

        p_matrix = matrix(p)
        for row in 1:nrows(p_matrix), column in 1:ncols(p_matrix)
            result[row, column] = p_matrix[row, column]
        end

        q_matrix = matrix(q)
        row_offset = ngens(domain(p))
        column_offset = ngens(codomain(p))
        for row in 1:nrows(q_matrix), column in 1:ncols(q_matrix)
            result[row_offset + row, column_offset + column] =
                q_matrix[row, column]
        end

        return result
    end

    for k in [1, 4], m in 0:2, n in 0:2

        context = MackeyContext(GAP.Globals.SymmetricGroup(k))
        zero_functor = constant_mackey_functor(context, free_module(ZZ, m))
        nonzero_functor = constant_mackey_functor(context, free_module(ZZ, n))

        zero_plus_nonzero = MackeyFunctors.direct_sum_mf(zero_functor, nonzero_functor)
        @test all(eachindex(context.subgroups)) do subgroup_index
            ngens(zero_plus_nonzero.values[subgroup_index]) == m + n
        end

        nonzero_plus_zero = MackeyFunctors.direct_sum_mf(nonzero_functor, zero_functor)
        @test all(eachindex(context.subgroups)) do subgroup_index
            ngens(nonzero_plus_zero.values[subgroup_index]) == m + n
        end

        for (sum_functor, left_functor, right_functor) in (
            (zero_plus_nonzero, zero_functor, nonzero_functor),
            (nonzero_plus_zero, nonzero_functor, zero_functor),
        )
            for (cover_index, (i, j)) in enumerate(context.covers)
                restriction_map = sum_functor.cover_restrictions[cover_index]
                @test matrix(restriction_map) == direct_sum_test_matrix(
                    left_functor.cover_restrictions[cover_index],
                    right_functor.cover_restrictions[cover_index],
                )

                transfer_map = sum_functor.cover_transfers[cover_index]
                @test matrix(transfer_map) == direct_sum_test_matrix(
                    left_functor.cover_transfers[cover_index],
                    right_functor.cover_transfers[cover_index],
                )
            end

            for generator_index in eachindex(context.generators)
                for subgroup_index in eachindex(context.subgroups)
                    conjugation_map =
                        sum_functor.generator_conjugations[generator_index, subgroup_index]
                    @test matrix(conjugation_map) == direct_sum_test_matrix(
                        left_functor.generator_conjugations[
                            generator_index,
                            subgroup_index,
                        ],
                        right_functor.generator_conjugations[
                            generator_index,
                            subgroup_index,
                        ],
                    )
                end
            end
        end
    end

    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    same_group_context = MackeyContext(context.group)
    @test MackeyFunctors.direct_sum_mf(
        constant_mackey_functor(context, ZZ),
        constant_mackey_functor(same_group_context, ZZ),
    ) isa MackeyFunctor

    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    Z2, = quo(F1, twoF1)

    z2_functor = constant_mackey_functor(context, Z2)
    z2_plus_z2 = MackeyFunctors.direct_sum_mf(z2_functor, z2_functor)
    @test all(eachindex(context.subgroups)) do subgroup_index
        AbstractAlgebra.invariant_factors(z2_plus_z2.values[subgroup_index]) ==
            BigInt[2, 2]
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

@testset "Zero Mackey functors" begin
    @test zero_mackey_functor(MackeyContext(GAP.Globals.CyclicGroup(2)), ZZ) isa MackeyFunctor
    @test coefficient_ring(
        zero_mackey_functor(MackeyContext(GAP.Globals.CyclicGroup(2)), GF(2)),
    ) == GF(2)
end
