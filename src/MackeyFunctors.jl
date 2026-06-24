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
    double_coset_representative_words,
    restriction,
    transfer

include("types/MackeyFunctor.jl")
export MackeyFunctor,
    coefficient_ring

include("Constructors.jl")
export constant_mackey_functor, burnside_mackey_functor

end
