module MackeyFunctors
using GAP
using AbstractAlgebra

# Import is needed so that our version doesn't clash with the one from AbstractAlgebra
import AbstractAlgebra: coefficient_ring

# We add some additional functionality needed from the abstract algebra world
include("abstract_algebra/ModuleHomomorphisms.jl")
include("abstract_algebra/HomModule.jl")
export HomModule, underlying_module, as_hom_module_element, as_homomorphism

include("group_theory/Words.jl")
export generator_relations

include("types/MackeyContext.jl")
export MackeyContext,
    double_coset_representative_data,
    double_coset_representative_words

include("GModules.jl")
export GModule, permutation_module

include("types/MackeyFunctor.jl")
export MackeyFunctor,
    coefficient_ring,
    conjugation,
    restriction,
    transfer,
    value

include("Show.jl")

include("Shift.jl")
export shift

include("Constructors.jl")
export constant_mackey_functor, burnside_mackey_functor,
    free_mackey_functor, fixedpoint_mackey_functor

include("types/Homomorphism.jl")
export MackeyFunctorHomomorphism, id_homomorphism

end
