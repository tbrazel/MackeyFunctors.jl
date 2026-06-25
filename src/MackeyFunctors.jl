module MackeyFunctors
using GAP
using AbstractAlgebra

# Importing this is needed so that our version doesn't clash with the one from AbstractAlgebra
import AbstractAlgebra: coefficient_ring

# We add some additional functionality needed from the abstract algebra world
include("abstract_algebra/ModuleHomomorphisms.jl")
include("abstract_algebra/DirectSums.jl")
include("abstract_algebra/HomModule.jl")
export HomModule, underlying_module, as_hom_module_element, as_homomorphism
include("abstract_algebra/TensorProduct.jl")
export TensorProduct,tensor_product,tensor_product_element

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
include("abstract_algebra/GModules.jl")
export GModule, permutation_module

# Various constructor methods
include("constructors/Constructors.jl")
export constant_mackey_functor, burnside_mackey_functor,
    free_mackey_functor, fixedpoint_mackey_functor,zero_mackey_functor

# Direct sum of Mackey functors and homomorphisms
include("DirectSum.jl")
export direct_sum_mf, direct_sum_homomorphism


end
