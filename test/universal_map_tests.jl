function _universal_map_with_element(
    M::MackeyFunctor,
    H_index::Int,
    x,
)
    burnside, shifted_burnside, basis_subgroup_indices =
        MackeyFunctors._shifted_burnside_mackey_functor(
            M.context,
            H_index,
            coefficient_ring(M),
        )
    universal = MackeyFunctors._universal_element(
        burnside,
        shifted_burnside,
        H_index,
        basis_subgroup_indices,
    )
    map = MackeyFunctors._universal_map(
        M,
        H_index,
        x,
        burnside,
        shifted_burnside,
        basis_subgroup_indices,
    )

    return map, universal
end

@testset "Universal elements and maps" begin
    C1 = GAP.Globals.CyclicGroup(1)
    c1_context = MackeyContext(C1)
    c1_constant = constant_mackey_functor(c1_context, ZZ)
    c1_x = ZZ(3) * gens(value(c1_constant, 1))[1]
    c1_A, c1_u = universal_element(c1_context, 1, ZZ)
    c1_public_map = universal_map(c1_constant, 1, c1_x)

    @test c1_A isa MackeyFunctor
    @test parent(c1_u) === value(c1_A, 1)
    @test c1_public_map isa MackeyFunctorHomomorphism
    @test c1_public_map.components[1](gens(value(c1_public_map.domain, 1))[1]) == c1_x

    C2 = GAP.Globals.CyclicGroup(2)
    c2_context = MackeyContext(C2)
    c2_trivial =
        findfirst(H -> Int(GAP.Globals.Size(H)) == 1, c2_context.subgroups)
    c2_constant = constant_mackey_functor(c2_context, ZZ)
    c2_x = ZZ(5) * gens(value(c2_constant, c2_trivial))[1]
    c2_map, c2_u =
        _universal_map_with_element(c2_constant, c2_trivial, c2_x)

    @test universal_map(c2_constant, c2_trivial, c2_x) isa MackeyFunctorHomomorphism
    @test c2_map.components[c2_trivial](c2_u) == c2_x

    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    Z2, = quo(F1, twoF1)
    c2_constant_Z2 = constant_mackey_functor(c2_context, Z2)
    c2_z2_x = gens(value(c2_constant_Z2, c2_trivial))[1]
    c2_z2_map, c2_z2_u =
        _universal_map_with_element(c2_constant_Z2, c2_trivial, c2_z2_x)

    @test c2_z2_map.components[c2_trivial](c2_z2_u) == c2_z2_x

    S3 = GAP.Globals.SymmetricGroup(3)
    s3_context = MackeyContext(S3)
    s3_order_two =
        findfirst(H -> Int(GAP.Globals.Size(H)) == 2, s3_context.subgroups)
    s3_burnside = burnside_mackey_functor(s3_context, ZZ)
    s3_x = gens(value(s3_burnside, s3_order_two))[1]
    s3_map, s3_u =
        _universal_map_with_element(s3_burnside, s3_order_two, s3_x)

    @test universal_map(s3_burnside, s3_order_two, s3_x) isa MackeyFunctorHomomorphism
    @test s3_map.components[s3_order_two](s3_u) == s3_x

    wrong_module = free_module(ZZ, 2)
    @test_throws ArgumentError universal_map(
        c2_constant,
        c2_trivial,
        gens(wrong_module)[1],
    )
end
