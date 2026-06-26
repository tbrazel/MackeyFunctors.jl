"""
    HomModule(M::FPModule, N::FPModule)

Represent ``Hom_R(M, N)`` as a finitely presented module over the common
base ring of `M` and `N`.

The underlying module is available as [`underlying_module`](@ref).  Elements
of that module can be converted to module homomorphisms with
[`as_homomorphism`](@ref), and module homomorphisms can be converted back with
[`as_hom_module_element`](@ref).
"""
struct HomModule{T <: RingElement}
    # The actual finitely presented module representing Hom_R(M, N).
    H::AbstractAlgebra.FPModule{T}

    # The source and target modules M and N. We keep these so that elements of
    # H can be turned back into AbstractAlgebra module homomorphisms M -> N.
    domain_module::AbstractAlgebra.FPModule{T}
    codomain_module::AbstractAlgebra.FPModule{T}

    # A homomorphism M -> N is determined by choosing one element of N for each
    # generator of M. Before imposing the relations of M, all such choices form
    # the ambient module N^(number of generators of M).
    ambient_module::AbstractAlgebra.FPModule{T}

    # For each relation of M, evaluating that relation on the chosen images
    # gives one element of N. These relation-values live in
    # N^(number of relations of M).
    relation_module::AbstractAlgebra.FPModule{T}

    # Sends a tuple of candidate generator images in ambient_module to the tuple
    # of relation-values in relation_module. The Hom module is its kernel.
    relation_map::Generic.ModuleHomomorphism{T}

    # The canonical inclusion H -> ambient_module returned by kernel.
    inclusion::Generic.ModuleHomomorphism{T}
end

# Return M^copies as an FPModule.
#
# HomModule relies on presentation data seen by kernel, so positive powers go
# through the package-local direct-sum workaround instead of
# AbstractAlgebra.direct_sum.  The zero-copy case is handled here because an
# empty direct sum needs the base ring of M to construct the zero-rank free
# module.
function _power_module(M::AbstractAlgebra.FPModule{T}, copies::Int) where T <: RingElement
    copies >= 0 || throw(ArgumentError("The number of copies must be nonnegative."))

    copies == 0 && return free_module(base_ring(M), 0)
    copies == 1 && return M

    return _direct_sum_module(AbstractAlgebra.FPModule{T}[M for _ in 1:copies])
end

function _relation_map(
    domain_module::AbstractAlgebra.FPModule{T},
    codomain_module::AbstractAlgebra.FPModule{T},
    ambient_module::AbstractAlgebra.FPModule{T},
    relation_module::AbstractAlgebra.FPModule{T},
) where T <: RingElement
    R = base_ring(domain_module)
    domain_generators = ngens(domain_module)
    codomain_generators = ngens(codomain_module)
    domain_relations = relations(domain_module)

    # AbstractAlgebra represents homomorphisms by matrices acting on row
    # vectors. This matrix sends a candidate tuple of generator images to the
    # tuple of evaluated relations.
    matrix_entries = zero_matrix(R, ngens(ambient_module), ngens(relation_module))
    for (relation_index, relation) in enumerate(domain_relations)
        # If relation a_1 e_1 + ... + a_m e_m = 0 holds in M, then a map
        # M -> N must satisfy a_1 f(e_1) + ... + a_m f(e_m) = 0 in N.
        # This matrix records all those linear conditions at once.
        for domain_generator in 1:domain_generators
            coefficient = relation[1, domain_generator]
            iszero(coefficient) && continue

            for codomain_generator in 1:codomain_generators
                # Coordinates are flattened by domain generator first:
                # (image of M-generator i, N-coordinate j).
                source_index = (domain_generator - 1)*codomain_generators + codomain_generator

                # The target has one copy of N for each relation of M.
                target_index = (relation_index - 1)*codomain_generators + codomain_generator
                matrix_entries[source_index, target_index] = coefficient
            end
        end
    end

    return ModuleHomomorphism(ambient_module, relation_module, matrix_entries)
end

function HomModule(
    domain_module::AbstractAlgebra.FPModule{T},
    codomain_module::AbstractAlgebra.FPModule{T},
) where T <: RingElement
    base_ring(domain_module) == base_ring(codomain_module) ||
        throw(ArgumentError("Modules must be defined over the same base ring."))

    # ambient_module is the space of all possible images of the generators of M.
    # For example, if M has three generators, this is N + N + N.
    ambient_module = _power_module(codomain_module, ngens(domain_module))

    # relation_module is where relation violations are recorded. One summand
    # of N is used for each defining relation of M.
    relation_module = _power_module(codomain_module, length(relations(domain_module)))

    # relation_map sends a candidate tuple of generator images to the tuple of
    # evaluated domain relations.
    relation_map = _relation_map(domain_module, codomain_module, ambient_module, relation_module)

    # Hom_R(M, N) consists exactly of the choices of images for the generators
    # of M that make every relation of M vanish in N.
    hom_module, inclusion = kernel(relation_map)

    return HomModule(
        hom_module,
        domain_module,
        codomain_module,
        ambient_module,
        relation_module,
        relation_map,
        inclusion,
    )
end

"""
    underlying_module(H::HomModule)

Return the finitely presented module representing ``Hom_R(M, N)``.
"""
underlying_module(H::HomModule) = H.H

# Make HomModule behave like the FPModule it wraps for the basic operations
# used elsewhere in the package.
AbstractAlgebra.base_ring(H::HomModule) = base_ring(underlying_module(H))
AbstractAlgebra.number_of_generators(H::HomModule) = ngens(underlying_module(H))
AbstractAlgebra.gens(H::HomModule) = gens(underlying_module(H))
AbstractAlgebra.gen(H::HomModule, i::Int) = gen(underlying_module(H), i)
AbstractAlgebra.relations(H::HomModule) = relations(underlying_module(H))

function _matrix_entries_for_hom_module_element(H::HomModule{T}, x::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
    parent(x) === underlying_module(H) ||
        throw(ArgumentError("Element does not belong to this Hom module."))

    # Elements of H live in the kernel object returned by AbstractAlgebra.
    # The inclusion translates such an element back to the ambient module
    # N^ngens(M).
    ambient_element = H.inclusion(x)
    domain_generators = ngens(H.domain_module)
    codomain_generators = ngens(H.codomain_module)

    ngens(H.ambient_module) == domain_generators*codomain_generators ||
        error("Internal Hom module presentation has incompatible generator count.")

    # The inclusion puts x back into the ambient module. Reading the blocks of
    # that element gives the matrix of the represented homomorphism.
    return T[
        ambient_element[(domain_generator - 1)*codomain_generators + codomain_generator]
        for domain_generator in 1:domain_generators
        for codomain_generator in 1:codomain_generators
    ]
end

"""
    as_homomorphism(H::HomModule, x::FPModuleElem)

Convert an element of the underlying Hom module to the represented module
homomorphism.
"""
function as_homomorphism(H::HomModule{T}, x::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
    R = base_ring(H)
    domain_generators = ngens(H.domain_module)
    codomain_generators = ngens(H.codomain_module)
    entries = _matrix_entries_for_hom_module_element(H, x)

    # AbstractAlgebra represents a hom M -> N by an ngens(M) x ngens(N) matrix:
    # row i is the coordinate vector for the image of the i-th generator of M.
    hom_matrix = matrix(R, domain_generators, codomain_generators, entries)

    return ModuleHomomorphism(H.domain_module, H.codomain_module, hom_matrix)
end

"""
    as_hom_module_element(H::HomModule, f::ModuleHomomorphism)

Convert a module homomorphism `f` into the corresponding element of the
underlying Hom module.
"""
function as_hom_module_element(
    H::HomModule{T},
    f::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism),
) where T <: RingElement
    domain(f) === H.domain_module ||
        throw(ArgumentError("The homomorphism has the wrong domain."))
    codomain(f) === H.codomain_module ||
        throw(ArgumentError("The homomorphism has the wrong codomain."))

    hom_matrix = matrix(f)
    domain_generators = ngens(H.domain_module)
    codomain_generators = ngens(H.codomain_module)
    size(hom_matrix) == (domain_generators, codomain_generators) ||
        throw(ArgumentError("The homomorphism matrix has incompatible dimensions."))

    # Interpret the matrix of f as a point of the ambient module N^ngens(M).
    # It represents an element of Hom(M, N) precisely when it lies in the
    # kernel submodule stored as H.H.
    entries = T[
        hom_matrix[domain_generator, codomain_generator]
        for domain_generator in 1:domain_generators
        for codomain_generator in 1:codomain_generators
    ]
    ambient_element = H.ambient_module(entries)

    # Lift through the kernel inclusion H -> ambient_module. Failure would mean
    # the supplied map matrix does not actually respect the relations of M.
    has_lift, hom_element = has_preimage_with_preimage(H.inclusion, ambient_element)
    has_lift ||
        throw(ArgumentError("The homomorphism does not send domain relations to zero."))

    return hom_element
end

function Base.show(io::IO, H::HomModule)
    print(io, "Hom module from ")
    print(io, H.domain_module)
    print(io, " to ")
    print(io, H.codomain_module)
end
