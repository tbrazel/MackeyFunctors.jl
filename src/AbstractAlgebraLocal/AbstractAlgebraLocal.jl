module AbstractAlgebraLocal
using AbstractAlgebra

include("fixes/DirectSums.jl")
export direct_sum, _direct_sum_module

include("fixes/ModuleHomomorphisms.jl")
export submodules_matrix

include("additional_functionality/Morphisms.jl")
export as_homomorphism,
    block_homomorphism,
    zero_homomorphism,
    identity_homomorphism,
    identity_isomorphism,
    is_identity_module_homomorphism

include("additional_functionality/HomModule.jl")
export Hom,
    HomModule,
    as_hom_module_element,
    postcomposition_map,
    precomposition_map,
    underlying_module

include("additional_functionality/TensorProduct.jl")
export TensorProduct, TensorProductElem, tensor_product, structure_map

end
