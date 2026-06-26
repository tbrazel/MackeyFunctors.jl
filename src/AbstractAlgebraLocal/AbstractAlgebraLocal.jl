module AbstractAlgebraLocal
using AbstractAlgebra

include("fixes/DirectSums.jl")
export _direct_sum, _direct_sum_module

include("fixes/ModuleHomomorphisms.jl")
export as_homomorphism,
    direct_sum_homomorphism,
    identity_homomorphism,
    identity_isomorphism,
    is_identity_module_homomorphism,
    is_isomorphism,
    is_zero_module_homomorphism,
    map_eq,
    submodules_matrix,
    zero_homomorphism

include("additional_functionality/HomModule.jl")
export HomModule, as_hom_module_element, underlying_module

include("additional_functionality/TensorProduct.jl")
export TensorProduct, tensor_product, tensor_product_element

end
