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