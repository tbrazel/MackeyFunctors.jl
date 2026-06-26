# Examples

This example constructs constant $\mathbb{Z}$-Mackey functors and Burnside Mackey functors for a few common groups:

```julia
using Pkg;Pkg.activate(".") # if you haven't already
using GAP, AbstractAlgebra, MackeyFunctors

# Computes the constant and Burnside C2 Mackey Functor
M = FreeModule(ZZ, 1)
C2 = GAP.Globals.CyclicGroup(2)
C2_context = MackeyContext(C2)
C2_Mackey_constant = constant_mackey_functor(C2_context, M)
C2_Burnside = burnside_mackey_functor(C2_context)

# Computes the constant and Burnside C4 Mackey Functor
C4 = GAP.Globals.CyclicGroup(4)
C4_context = MackeyContext(C4)
C4_Mackey_constant = constant_mackey_functor(C4_context, M)
C4_Burnside = burnside_mackey_functor(C4_context)

# Computes the constant and Burnside S3 Mackey Functor
S3 = GAP.Globals.SymmetricGroup(3)
S3_context = MackeyContext(S3)
S3_Mackey_constant = constant_mackey_functor(S3_context, M)
S3_Burnside = burnside_mackey_functor(S3_context)

# Constructs a MackeyContext for S5
S5 = GAP.Globals.SymmetricGroup(5)
S5_context = MackeyContext(S5)

```