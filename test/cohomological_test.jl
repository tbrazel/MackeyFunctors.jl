@testset "Cohomological Burnside Boolean" begin

    C1 = GAP.Globals.CyclicGroup(1)
    C1_context = MackeyContext(C1)
    burnside_C1 = burnside_mackey_functor(C1_context)
    @test is_cohomological(burnside_C1) == true

    C3 = GAP.Globals.CyclicGroup(3)
    C3_context = MackeyContext(C3)
    burnside_C3 = burnside_mackey_functor(C3_context)
    @test is_cohomological(burnside_C3) == false

    V4 = GAP.Globals.DirectProduct(GAP.Globals.CyclicGroup(2), GAP.Globals.CyclicGroup(2))
    V4_context = MackeyContext(V4)
    burnside_V4 = burnside_mackey_functor(V4_context)
    @test is_cohomological(burnside_V4) == false

end

@testset "Cohomological Boolean C2" begin
    G = GAP.Globals.CyclicGroup(2)
    ctx = MackeyContext(G)

    # Zero functor 
    Z = zero_mackey_functor(ctx)
    @test is_cohomological(Z) == true

     # Constant functor 
    constZ = constant_mackey_functor(ctx, ZZ)
    @test is_cohomological(constZ) == true

    # Direct sum closure
    @test is_cohomological(direct_sum(constZ, Z)[1]) == true

end