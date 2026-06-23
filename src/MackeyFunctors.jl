module MackeyFunctors
using GAP
using AbstractAlgebra

include("MiscFuncs.jl")

include("module_methods/Modules.jl")

include("group_theory/AbelianGroupHom.jl")
export defines_fg_abelian_map, are_equal_abelian_group_homomorphisms # don't need to export these long term

include("group_theory/Words.jl")
export generator_relations

include("group_theory/HelperFunctions.jl")

include("types/MackeyContext.jl")
export MackeyContext, double_coset_representative_data, double_coset_representative_words

include("types/MackeyFunctor.jl")
export MackeyFunctor

include("group_theory/SubgroupLattice.jl")
export all_supergroups_from_lattice

include("visualizers/MackeyFunctorVisualizerData.jl")
export visualizer_data, visualizer_json

end
