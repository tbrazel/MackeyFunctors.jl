@testset "Internal Hom Mackey functors" begin
    C1 = GAP.Globals.CyclicGroup(1)
    c1_context = MackeyContext(C1)
    c1_M = constant_mackey_functor(c1_context, ZZ)
    c1_N = constant_mackey_functor(c1_context, ZZ)
    c1_internal_hom = InternalHom(c1_M, c1_N)
    c1_external_hom = Hom(c1_M, c1_N)

    @test c1_internal_hom isa MackeyFunctor
    @test MackeyFunctorHomMackeyFunctor(c1_M, c1_N) isa MackeyFunctor
    @test ngens(value(c1_internal_hom, 1)) ==
        ngens(underlying_module(c1_external_hom))
    @test isempty(c1_internal_hom.cover_restrictions)
    @test isempty(c1_internal_hom.cover_transfers)

    C4 = GAP.Globals.CyclicGroup(4)
    c4_context = MackeyContext(C4)
    c4_M = burnside_mackey_functor(c4_context)
    c4_N = constant_mackey_functor(c4_context, ZZ)
    c4_internal_hom = InternalHom(c4_M, c4_N)
    c4_hom_modules = MackeyFunctorHomModule[
        Hom(shift(c4_M, H_index), c4_N)
        for H_index in eachindex(c4_context.subgroups)
    ]

    @test c4_internal_hom isa MackeyFunctor
    @test all(eachindex(c4_context.subgroups)) do H_index
        ngens(value(c4_internal_hom, H_index)) ==
            ngens(underlying_module(c4_hom_modules[H_index]))
    end

    for (cover_index, (H_index, K_index)) in enumerate(c4_context.covers)
        expected_restriction = precomposition_map(
            c4_hom_modules[K_index],
            c4_hom_modules[H_index],
            shift_transfer(c4_M, H_index, K_index),
        )
        expected_transfer = precomposition_map(
            c4_hom_modules[H_index],
            c4_hom_modules[K_index],
            shift_restriction(c4_M, H_index, K_index),
        )

        @test matrix(c4_internal_hom.cover_restrictions[cover_index]) ==
            matrix(expected_restriction)
        @test matrix(c4_internal_hom.cover_transfers[cover_index]) ==
            matrix(expected_transfer)
    end

    S3 = GAP.Globals.SymmetricGroup(3)
    s3_context = MackeyContext(S3)
    s3_M = constant_mackey_functor(s3_context, ZZ)
    s3_N = constant_mackey_functor(s3_context, ZZ)
    s3_internal_hom = InternalHom(s3_M, s3_N)
    s3_hom_modules = MackeyFunctorHomModule[
        Hom(shift(s3_M, H_index), s3_N)
        for H_index in eachindex(s3_context.subgroups)
    ]

    for generator_index in eachindex(s3_context.generators), H_index in eachindex(s3_context.subgroups)
        g = s3_context.generators[generator_index]
        target_H_index =
            s3_context.generator_left_conjugation_matrix[generator_index, H_index]
        expected_conjugation = precomposition_map(
            s3_hom_modules[H_index],
            s3_hom_modules[target_H_index],
            shift_conjugation(s3_M, g^-1, target_H_index),
        )

        @test matrix(s3_internal_hom.generator_conjugations[generator_index, H_index]) ==
            matrix(expected_conjugation)
    end

    c2_context = MackeyContext(GAP.Globals.CyclicGroup(2))
    c3_context = MackeyContext(GAP.Globals.CyclicGroup(3))
    @test_throws ArgumentError InternalHom(
        constant_mackey_functor(c2_context, ZZ),
        constant_mackey_functor(c3_context, ZZ),
    )
    @test_throws ArgumentError InternalHom(
        constant_mackey_functor(c2_context, ZZ),
        constant_mackey_functor(c2_context, GF(2)),
    )
end
