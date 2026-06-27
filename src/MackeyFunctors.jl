module MackeyFunctors
using GAP
using AbstractAlgebra

include("AbstractAlgebraLocal/AbstractAlgebraLocal.jl")
using .AbstractAlgebraLocal
using .AbstractAlgebraLocal: ModuleHomomorphism, ModuleIsomorphism, is_invertible

# Importing this is needed so that our version doesn't clash with the one from AbstractAlgebra
import AbstractAlgebra: coefficient_ring

# Some methods for manipulating words in generators of a group
include("group_theory/Words.jl")
export generator_relations

# Defines our MackeyContext type
include("types/MackeyContext.jl")
export MackeyContext

# MackeyFunctor type
include("types/MackeyFunctor.jl")
export MackeyFunctor,
    coefficient_ring,
    conjugation,
    restriction,
    transfer,
    value

# MackeyFunctorHomomorphism type
include("types/Homomorphism.jl")
export MackeyFunctorHomomorphism, id_homomorphism

# Nice printing for various new types
include("Show.jl")

# Provides the "shift" operation which helps us construct new Mackey functors out of old ones
include("constructors/Shift.jl")
export shift

# Basic theory of G-modules
include("types/GModules.jl")
export GModule, permutation_module

# Various constructor methods
include("constructors/Constructors.jl")
export constant_mackey_functor, burnside_mackey_functor,
    free_mackey_functor, fixedpoint_mackey_functor,zero_mackey_functor

# Direct sum of Mackey functors and homomorphisms
include("constructors/DirectSum.jl")
export direct_sum, direct_sum_homomorphism

# Cohomological Mackey functors
include("Cohomological.jl")
export is_cohomological

end
