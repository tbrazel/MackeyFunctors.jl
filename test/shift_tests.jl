using MackeyFunctors.AbstractAlgebraLocal: identity_homomorphism

@testset "Shifting identity homomorphisms on Mackey functors" begin
    # For S3
    G = GAP.Globals.SymmetricGroup(3)
    ctx = MackeyContext(G)

    # Make the Burnside Mackey functor for S3 and its identity homomorphism
    A = burnside_mackey_functor(ctx)
    idA = id_homomorphism(A)

    # Assert that the components of the shift of the identity map are still identity maps
    for grp_index in eachindex(ctx.subgroups)
        idA_H = shift(idA,grp_index)
        @test idA_H isa MackeyFunctorHomomorphism
        for k in eachindex(ctx.subgroups)
            @test idA_H.components[k] == identity_homomorphism(idA_H.domain.values[k])
        end
    end

end

@testset "Natural maps between shifted Mackey functors" begin
    C4 = GAP.Globals.CyclicGroup(4)
    c4_context = MackeyContext(C4)
    c4_burnside = burnside_mackey_functor(c4_context)

    c4_trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, c4_context.subgroups)
    c4_order_two = findfirst(H -> Int(GAP.Globals.Size(H)) == 2, c4_context.subgroups)
    c4_whole = findfirst(
        H -> Int(GAP.Globals.Size(H)) == Int(GAP.Globals.Size(C4)),
        c4_context.subgroups,
    )

    identity_transfer = shift_transfer(c4_burnside, c4_trivial, c4_trivial)
    identity_restriction = shift_restriction(c4_burnside, c4_trivial, c4_trivial)
    @test identity_transfer isa MackeyFunctorHomomorphism
    @test identity_restriction isa MackeyFunctorHomomorphism
    @test all(eachindex(c4_context.subgroups)) do K_index
        identity_transfer.components[K_index] ==
            identity_homomorphism(identity_transfer.domain.values[K_index])
    end
    @test all(eachindex(c4_context.subgroups)) do K_index
        identity_restriction.components[K_index] ==
            identity_homomorphism(identity_restriction.domain.values[K_index])
    end

    trivial_to_order_two = shift_transfer(c4_burnside, c4_trivial, c4_order_two)
    order_two_to_whole = shift_transfer(c4_burnside, c4_order_two, c4_whole)
    trivial_to_whole = shift_transfer(c4_burnside, c4_trivial, c4_whole)
    @test trivial_to_order_two isa MackeyFunctorHomomorphism
    @test order_two_to_whole isa MackeyFunctorHomomorphism
    @test all(eachindex(c4_context.subgroups)) do K_index
        trivial_to_order_two.components[K_index] *
        order_two_to_whole.components[K_index] ==
            trivial_to_whole.components[K_index]
    end

    order_two_to_trivial = shift_restriction(c4_burnside, c4_trivial, c4_order_two)
    whole_to_order_two = shift_restriction(c4_burnside, c4_order_two, c4_whole)
    whole_to_trivial = shift_restriction(c4_burnside, c4_trivial, c4_whole)
    @test order_two_to_trivial isa MackeyFunctorHomomorphism
    @test whole_to_order_two isa MackeyFunctorHomomorphism
    @test all(eachindex(c4_context.subgroups)) do K_index
        whole_to_order_two.components[K_index] *
        order_two_to_trivial.components[K_index] ==
            whole_to_trivial.components[K_index]
    end

    @test_throws ArgumentError shift_transfer(c4_burnside, c4_whole, c4_trivial)
    @test_throws ArgumentError shift_restriction(c4_burnside, c4_whole, c4_trivial)

    S3 = GAP.Globals.SymmetricGroup(3)
    s3_context = MackeyContext(S3)
    s3_burnside = burnside_mackey_functor(s3_context)
    s3_order_two = findfirst(H -> Int(GAP.Globals.Size(H)) == 2, s3_context.subgroups)
    s3_H = s3_context.subgroups[s3_order_two]

    non_normalizing_element = nothing
    for candidate in GAP.Globals.Elements(S3)
        if s3_H^(candidate^-1) != s3_H
            non_normalizing_element = candidate
            break
        end
    end
    @test non_normalizing_element !== nothing

    identity_element = GAP.Globals.One(S3)
    identity_conjugation = shift_conjugation(s3_burnside, identity_element, s3_order_two)
    @test identity_conjugation isa MackeyFunctorHomomorphism
    @test all(eachindex(s3_context.subgroups)) do K_index
        identity_conjugation.components[K_index] ==
            identity_homomorphism(identity_conjugation.domain.values[K_index])
    end

    g = non_normalizing_element
    target_order_two = MackeyFunctors.subgroup_index(s3_context, s3_H^(g^-1))
    conjugation_by_g = shift_conjugation(s3_burnside, g, s3_order_two)
    conjugation_by_g_inverse = shift_conjugation(s3_burnside, g^-1, target_order_two)
    @test conjugation_by_g isa MackeyFunctorHomomorphism
    @test conjugation_by_g.codomain === shift(s3_burnside, target_order_two)
    @test conjugation_by_g_inverse isa MackeyFunctorHomomorphism
    @test conjugation_by_g_inverse.codomain === shift(s3_burnside, s3_order_two)
    @test all(eachindex(s3_context.subgroups)) do K_index
        conjugation_by_g.components[K_index] *
        conjugation_by_g_inverse.components[K_index] ==
            identity_homomorphism(conjugation_by_g.domain.values[K_index])
    end
end
