module MackeyFunctors
using GAP
using AbstractAlgebra

# Import is needed so that our version doesn't clash with the one from AbstractAlgebra
import AbstractAlgebra: coefficient_ring

include("module_methods/Modules.jl")

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
