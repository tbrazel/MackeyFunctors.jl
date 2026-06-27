module AbstractAlgebraLocal
using AbstractAlgebra

include("fixes/DirectSums.jl")
export _direct_sum, _direct_sum_module

include("fixes/ModuleHomomorphisms.jl")
export direct_sum_homomorphism,
    submodules_matrix

include("additional_functionality/Morphisms.jl")
export as_homomorphism,
    zero_homomorphism,
    identity_homomorphism,
    identity_isomorphism,
    is_identity_module_homomorphism

include("additional_functionality/HomModule.jl")
export HomModule, as_hom_module_element, underlying_module

include("additional_functionality/TensorProduct.jl")
export TensorProduct, tensor_product, tensor_product_element

end
