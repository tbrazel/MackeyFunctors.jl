using Pkg; Pkg.activate("Documents/GitHub/MackeyFunctors.jl")
using GAP, AbstractAlgebra, MackeyFunctors

M = FreeModule(ZZ, 1)
C2 = GAP.Globals.CyclicGroup(2)
C2_context = MackeyContext(C2)
C2_Mackey_constant = constant_mackey_functor(C2_context, M)
C2_Burnside = burnside_mackey_functor(C2_context)


C4 = GAP.Globals.CyclicGroup(4)
C4_context = MackeyContext(C4)
C4_Mackey_constant = constant_mackey_functor(C4_context, M)
C4_Burnside = burnside_mackey_functor(C4_context)


S3 = GAP.Globals.SymmetricGroup(3)
S3_context = MackeyContext(S3)
S3_Mackey_constant = constant_mackey_functor(S3_context, M)
S3_Burnside = burnside_mackey_functor(S3_context)


S5 = GAP.Globals.SymmetricGroup(5)
S5_context = MackeyContext(S5)





