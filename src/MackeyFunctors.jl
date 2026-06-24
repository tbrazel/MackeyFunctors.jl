module MackeyFunctors
using GAP
using AbstractAlgebra

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
export MackeyFunctor

include("Constructors.jl")
export constant_mackey_functor

end
