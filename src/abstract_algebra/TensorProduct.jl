"""
    TensorProduct(M::FPModule, N::FPModule)

Represent ``M \\otimes_R N`` as a finitely presented module over the common
base ring of `M` and `N`.

The underlying module is available as [`underlying_module`](@ref), and pure
tensors can be constructed with [`tensor_product_element`](@ref).
"""
struct TensorProduct{T <: RingElement}
    # The quotient FPModule representing M tensor_R N.
    tensor_module::AbstractAlgebra.FPModule{T}

    # The original factors are kept so that pure tensors and induced maps can
    # validate their inputs and recover the intended presentation.
    left_factor::AbstractAlgebra.FPModule{T}
    right_factor::AbstractAlgebra.FPModule{T}

    # Before imposing relations, the tensor product has one free generator
    # e_i tensor f_j for every pair of generators of the two factors.
    free_tensor_module::AbstractAlgebra.FPModule{T}

    # The quotient map from that free module to tensor_module.  We keep it
    # because AbstractAlgebra may simplify the quotient presentation, so the
    # generators of tensor_module need not be the original pair-generators.
    projection_from_free::Generic.ModuleHomomorphism{T}
end

function _zero_tensor_relation_row(R::Ring, n::Int)
    return [zero(R) for _ in 1:n]
end

function _is_zero_tensor_relation_row(row)
    return all(iszero, row)
end

function _tensor_generator_index(right_generator_count::Int, left_generator::Int, right_generator::Int)
    # Pair-generators are stored in row-major order:
    # (e_1 tensor f_1), ..., (e_1 tensor f_n),
    # (e_2 tensor f_1), ...
    return (left_generator - 1)*right_generator_count + right_generator
end

function _tensor_product_module(
    left_factor::AbstractAlgebra.FPModule{T},
    right_factor::AbstractAlgebra.FPModule{T},
) where T <: RingElement
    base_ring(left_factor) == base_ring(right_factor) ||
        throw(ArgumentError("Modules must be defined over the same base ring."))

    R = base_ring(left_factor)
    left_generators = ngens(left_factor)
    right_generators = ngens(right_factor)
    tensor_generators = left_generators*right_generators

    # This is the ambient free module on formal pair-generators.  The two
    # relation loops below quotient out exactly the relations forced by the
    # presentations of the two input modules.
    free_tensor_module = free_module(R, tensor_generators)

    relation_generators = elem_type(free_tensor_module)[]

    # If sum_i a_i e_i = 0 is a relation in M, impose
    # sum_i a_i (e_i tensor f_j) = 0 for every generator f_j of N.
    for relation in relations(left_factor), right_generator in 1:right_generators
        row = _zero_tensor_relation_row(R, tensor_generators)
        for left_generator in 1:left_generators
            row[
                _tensor_generator_index(
                    right_generators,
                    left_generator,
                    right_generator,
                )
            ] = relation[1, left_generator]
        end
        _is_zero_tensor_relation_row(row) ||
            push!(relation_generators, free_tensor_module(row))
    end

    # If sum_j b_j f_j = 0 is a relation in N, impose
    # sum_j b_j (e_i tensor f_j) = 0 for every generator e_i of M.
    for left_generator in 1:left_generators, relation in relations(right_factor)
        row = _zero_tensor_relation_row(R, tensor_generators)
        for right_generator in 1:right_generators
            row[
                _tensor_generator_index(
                    right_generators,
                    left_generator,
                    right_generator,
                )
            ] = relation[1, right_generator]
        end
        _is_zero_tensor_relation_row(row) ||
            push!(relation_generators, free_tensor_module(row))
    end

    if isempty(relation_generators)
        # Avoid taking a quotient by the zero relation submodule.  Returning the
        # free module with its identity map keeps the same interface while
        # preserving a cleaner presentation for free tensor products.
        return free_tensor_module, free_tensor_module, identity_homomorphism(free_tensor_module)
    end

    relation_submodule, = sub(free_tensor_module, relation_generators)
    tensor_module, projection_from_free = quo(free_tensor_module, relation_submodule)

    return tensor_module, free_tensor_module, projection_from_free
end

function TensorProduct(
    left_factor::AbstractAlgebra.FPModule{T},
    right_factor::AbstractAlgebra.FPModule{T},
) where T <: RingElement
    tensor_module, free_tensor_module, projection_from_free =
        _tensor_product_module(left_factor, right_factor)

    return TensorProduct(
        tensor_module,
        left_factor,
        right_factor,
        free_tensor_module,
        projection_from_free,
    )
end

function TensorProduct(
    left_factor::AbstractAlgebra.FPModule,
    right_factor::AbstractAlgebra.FPModule,
)
    # This fallback gives a clearer error when two modules are both FPModules
    # but do not have the same coefficient element type, so the parametric
    # constructor above cannot apply.
    base_ring(left_factor) == base_ring(right_factor) ||
        throw(ArgumentError("Modules must be defined over the same base ring."))

    throw(ArgumentError("Modules must have compatible coefficient element types."))
end

"""
    tensor_product(M::FPModule, N::FPModule)

Construct the tensor product wrapper representing ``M \\otimes_R N``.
"""
tensor_product(
    left_factor::AbstractAlgebra.FPModule,
    right_factor::AbstractAlgebra.FPModule,
) = TensorProduct(left_factor, right_factor)

"""
    underlying_module(T::TensorProduct)

Return the finitely presented module representing the tensor product.
"""
underlying_module(T::TensorProduct) = T.tensor_module

# Make TensorProduct behave like the FPModule it wraps for the basic operations
# used elsewhere in the package.
AbstractAlgebra.base_ring(T::TensorProduct) = base_ring(underlying_module(T))
AbstractAlgebra.number_of_generators(T::TensorProduct) = ngens(underlying_module(T))
AbstractAlgebra.gens(T::TensorProduct) = gens(underlying_module(T))
AbstractAlgebra.gen(T::TensorProduct, i::Int) = gen(underlying_module(T), i)
AbstractAlgebra.relations(T::TensorProduct) = relations(underlying_module(T))

"""
    tensor_product_element(T::TensorProduct, m::FPModuleElem, n::FPModuleElem)

Return the pure tensor `m tensor n` as an element of `underlying_module(T)`.
"""
function tensor_product_element(
    T::TensorProduct{S},
    left_element::AbstractAlgebra.FPModuleElem{S},
    right_element::AbstractAlgebra.FPModuleElem{S},
) where S <: RingElement
    parent(left_element) === T.left_factor ||
        throw(ArgumentError("Left element does not belong to the left tensor factor."))
    parent(right_element) === T.right_factor ||
        throw(ArgumentError("Right element does not belong to the right tensor factor."))

    R = base_ring(T)
    right_generators = ngens(T.right_factor)
    entries = [zero(R) for _ in 1:ngens(T.free_tensor_module)]

    # Expand m and n in the chosen generators and use bilinearity:
    # (sum_i a_i e_i) tensor (sum_j b_j f_j)
    #   = sum_{i,j} a_i*b_j (e_i tensor f_j).
    for left_generator in 1:ngens(T.left_factor), right_generator in 1:right_generators
        index = _tensor_generator_index(
            right_generators,
            left_generator,
            right_generator,
        )
        entries[index] += left_element[left_generator]*right_element[right_generator]
    end

    return T.projection_from_free(T.free_tensor_module(entries))
end

function _tensor_homomorphism_from_generator_images(
    domain_module::AbstractAlgebra.FPModule{T},
    codomain_module::AbstractAlgebra.FPModule{T},
    images::Vector{<:AbstractAlgebra.FPModuleElem{T}},
) where T <: RingElement
    if ngens(domain_module) == 0
        return ModuleHomomorphism(
            domain_module,
            codomain_module,
            zero_matrix(base_ring(domain_module), 0, ngens(codomain_module)),
        )
    end

    return ModuleHomomorphism(domain_module, codomain_module, images)
end

"""
    tensor_product(f::ModuleHomomorphism, g::ModuleHomomorphism)

Return the induced homomorphism
``domain(f) \\otimes domain(g) -> codomain(f) \\otimes codomain(g)``.
"""
function tensor_product(
    left_map::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism),
    right_map::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism),
)
    source_tensor = TensorProduct(domain(left_map), domain(right_map))
    target_tensor = TensorProduct(codomain(left_map), codomain(right_map))
    source_module = underlying_module(source_tensor)
    target_module = underlying_module(target_tensor)

    # Given generators e_i for domain(f) and f_j for domain(g), this stores f(e_i) otimes g(f_j) as an element of the target tensor product
    pair_images = elem_type(target_module)[
        tensor_product_element(
            target_tensor,
            left_map(gen(domain(left_map), left_generator)),
            right_map(gen(domain(right_map), right_generator)),
        )
        for left_generator in 1:ngens(domain(left_map))
        for right_generator in 1:ngens(domain(right_map))
    ]

    generator_images = elem_type(target_module)[]
    for source_generator in gens(source_module)
        # Lift each actual generator of the source quotient back to the free
        # pair-generator presentation, then evaluate the corresponding linear
        # combination of pair_images.
        lift = preimage(source_tensor.projection_from_free, source_generator)
        image = zero(target_module)
        for index in 1:ngens(source_tensor.free_tensor_module)
            coefficient = lift[index]
            iszero(coefficient) && continue
            image += coefficient*pair_images[index]
        end
        push!(generator_images, image)
    end

    return _tensor_homomorphism_from_generator_images(
        source_module,
        target_module,
        generator_images,
    )
end

function Base.show(io::IO, T::TensorProduct)
    print(io, "Tensor product of ")
    print(io, T.left_factor)
    print(io, " and ")
    print(io, T.right_factor)
end
