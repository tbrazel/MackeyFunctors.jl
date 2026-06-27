@testset "Cohomological Boolean" begin

    C1 = GAP.Globals.CyclicGroup(1)
    C1_context = MackeyContext(C1)
    burnside_C1 = burnside_mackey_functor(C1_context)
    @test is_cohomological(burnside_C1) == true

end
