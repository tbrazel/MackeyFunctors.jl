using Test
using GAP
using AbstractAlgebra
using MackeyFunctors

using MackeyFunctors: ModuleHomomorphism, ModuleIsomorphism

include("shift_tests.jl")
include("abstract_algebra_tests.jl")
include("mackey_context_tests.jl")
include("mackey_functor_tests.jl")
include("homomorphism_tests.jl")
include("internal_hom_tests.jl")
include("universal_map_tests.jl")
include("cohomological_test.jl")
