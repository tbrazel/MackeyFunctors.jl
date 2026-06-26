"""
    MackeyFunctorHomModule(M::MackeyFunctor, N::MackeyFunctor)
    HomModule(M::MackeyFunctor, N::MackeyFunctor)

Represent the module of Mackey functor homomorphisms ``M -> N`` as a
finitely presented module over the common coefficient ring.

The underlying module is available as [`underlying_module`](@ref).  Elements
of that module can be converted to [`MackeyFunctorHomomorphism`](@ref)s with
[`as_homomorphism`](@ref), and Mackey functor homomorphisms can be converted
back with [`as_hom_module_element`](@ref).
"""
struct MackeyFunctorHomModule{T <: RingElement}
    # The actual finitely presented module representing Hom(M, N).
    H::AbstractAlgebra.FPModule{T}

    domain_mf::MackeyFunctor
    codomain_mf::MackeyFunctor

    # The unrestricted product of module Homs, one summand for each subgroup.
    component_hom_modules::Vector{HomModule{T}}
    ambient_module::AbstractAlgebra.FPModule{T}
    ambient_injections::Vector{Generic.ModuleHomomorphism{T}}
    ambient_projections::Vector{Generic.ModuleHomomorphism{T}}

    # The codomain of the compatibility map is a direct sum of module Homs,
    # one summand for every restriction, transfer, and conjugation equation.
    compatibility_hom_modules::Vector{HomModule{T}}
    compatibility_module::AbstractAlgebra.FPModule{T}
    compatibility_map::Generic.ModuleHomomorphism{T}

    # The canonical inclusion H -> ambient_module returned by kernel.
    inclusion::Generic.ModuleHomomorphism{T}
end

function _check_hom_module_map_rings(source::HomModule{T}, target::HomModule{T}) where T <: RingElement
    base_ring(source) == base_ring(target) ||
        throw(ArgumentError("Hom modules must be defined over the same base ring."))
    return nothing
end

function _hom_module_precomposition_map(
    source::HomModule{T},
    target::HomModule{T},
    p,
) where T <: RingElement
    _check_hom_module_map_rings(source, target)
    p = _as_homomorphism(p)

    domain(p) === target.domain_module ||
        throw(ArgumentError("The precomposition map has the wrong domain."))
    codomain(p) === source.domain_module ||
        throw(ArgumentError("The precomposition map has the wrong codomain."))
    source.codomain_module === target.codomain_module ||
        throw(ArgumentError("Precomposition target Hom module has the wrong codomain."))

    images = elem_type(underlying_module(target))[
        as_hom_module_element(target, p * as_homomorphism(source, x))
        for x in gens(source)
    ]

    return _homomorphism_from_generator_images(
        underlying_module(source),
        underlying_module(target),
        images,
    )
end

function _hom_module_postcomposition_map(
    source::HomModule{T},
    target::HomModule{T},
    q,
) where T <: RingElement
    _check_hom_module_map_rings(source, target)
    q = _as_homomorphism(q)

    domain(q) === source.codomain_module ||
        throw(ArgumentError("The postcomposition map has the wrong domain."))
    codomain(q) === target.codomain_module ||
        throw(ArgumentError("The postcomposition map has the wrong codomain."))
    source.domain_module === target.domain_module ||
        throw(ArgumentError("Postcomposition target Hom module has the wrong domain."))

    images = elem_type(underlying_module(target))[
        as_hom_module_element(target, as_homomorphism(source, x) * q)
        for x in gens(source)
    ]

    return _homomorphism_from_generator_images(
        underlying_module(source),
        underlying_module(target),
        images,
    )
end

function _push_compatibility_equation!(
    compatibility_hom_modules::Vector{HomModule{T}},
    left_component_indices::Vector{Int},
    left_maps::Vector{Generic.ModuleHomomorphism{T}},
    right_component_indices::Vector{Int},
    right_maps::Vector{Generic.ModuleHomomorphism{T}},
    equation_hom_module::HomModule{T},
    left_component_index::Int,
    left_map::Generic.ModuleHomomorphism{T},
    right_component_index::Int,
    right_map::Generic.ModuleHomomorphism{T},
) where T <: RingElement
    push!(compatibility_hom_modules, equation_hom_module)
    push!(left_component_indices, left_component_index)
    push!(left_maps, left_map)
    push!(right_component_indices, right_component_index)
    push!(right_maps, right_map)
    return nothing
end

function _mackey_functor_hom_compatibility_data(
    domain_mf::MackeyFunctor,
    codomain_mf::MackeyFunctor,
    component_hom_modules::Vector{HomModule{T}},
) where T <: RingElement
    context = domain_mf.context

    compatibility_hom_modules = HomModule{T}[]
    left_component_indices = Int[]
    right_component_indices = Int[]
    left_maps = Generic.ModuleHomomorphism{T}[]
    right_maps = Generic.ModuleHomomorphism{T}[]

    for (cover_index, (i, j)) in enumerate(context.covers)
        transfer_hom = HomModule(domain_mf.values[i], codomain_mf.values[j])
        transfer_left_map = _hom_module_postcomposition_map(
            component_hom_modules[i],
            transfer_hom,
            codomain_mf.cover_transfers[cover_index],
        )
        transfer_right_map = _hom_module_precomposition_map(
            component_hom_modules[j],
            transfer_hom,
            domain_mf.cover_transfers[cover_index],
        )
        _push_compatibility_equation!(
            compatibility_hom_modules,
            left_component_indices,
            left_maps,
            right_component_indices,
            right_maps,
            transfer_hom,
            i,
            transfer_left_map,
            j,
            transfer_right_map,
        )

        restriction_hom = HomModule(domain_mf.values[j], codomain_mf.values[i])
        restriction_left_map = _hom_module_precomposition_map(
            component_hom_modules[i],
            restriction_hom,
            domain_mf.cover_restrictions[cover_index],
        )
        restriction_right_map = _hom_module_postcomposition_map(
            component_hom_modules[j],
            restriction_hom,
            codomain_mf.cover_restrictions[cover_index],
        )
        _push_compatibility_equation!(
            compatibility_hom_modules,
            left_component_indices,
            left_maps,
            right_component_indices,
            right_maps,
            restriction_hom,
            i,
            restriction_left_map,
            j,
            restriction_right_map,
        )
    end

    for generator_index in eachindex(context.generators), subgroup_index in eachindex(context.subgroups)
        target_subgroup_index =
            context.generator_left_conjugation_matrix[generator_index, subgroup_index]

        conjugation_hom = HomModule(
            domain_mf.values[subgroup_index],
            codomain_mf.values[target_subgroup_index],
        )
        conjugation_left_map = _hom_module_precomposition_map(
            component_hom_modules[target_subgroup_index],
            conjugation_hom,
            domain_mf.generator_conjugations[generator_index, subgroup_index],
        )
        conjugation_right_map = _hom_module_postcomposition_map(
            component_hom_modules[subgroup_index],
            conjugation_hom,
            codomain_mf.generator_conjugations[generator_index, subgroup_index],
        )
        _push_compatibility_equation!(
            compatibility_hom_modules,
            left_component_indices,
            left_maps,
            right_component_indices,
            right_maps,
            conjugation_hom,
            target_subgroup_index,
            conjugation_left_map,
            subgroup_index,
            conjugation_right_map,
        )
    end

    return (
        compatibility_hom_modules,
        left_component_indices,
        left_maps,
        right_component_indices,
        right_maps,
    )
end

function _compatibility_module_and_injections(
    R::Ring,
    compatibility_hom_modules::Vector{HomModule{T}},
) where T <: RingElement
    if isempty(compatibility_hom_modules)
        return (
            free_module(R, 0),
            Generic.ModuleHomomorphism{T}[],
        )
    end

    compatibility_modules = AbstractAlgebra.FPModule{T}[
        underlying_module(H)
        for H in compatibility_hom_modules
    ]
    compatibility_module, compatibility_injections, _ =
        _direct_sum(compatibility_modules)

    return compatibility_module, compatibility_injections
end

function _mackey_functor_hom_compatibility_map(
    ambient_module::AbstractAlgebra.FPModule{T},
    ambient_projections::Vector{Generic.ModuleHomomorphism{T}},
    compatibility_module::AbstractAlgebra.FPModule{T},
    compatibility_injections::Vector{Generic.ModuleHomomorphism{T}},
    left_component_indices::Vector{Int},
    left_maps::Vector{Generic.ModuleHomomorphism{T}},
    right_component_indices::Vector{Int},
    right_maps::Vector{Generic.ModuleHomomorphism{T}},
) where T <: RingElement
    images = elem_type(compatibility_module)[]

    for ambient_generator in gens(ambient_module)
        compatibility_value = zero(compatibility_module)

        for equation_index in eachindex(left_maps)
            left_value = left_maps[equation_index](
                ambient_projections[left_component_indices[equation_index]](
                    ambient_generator,
                ),
            )
            right_value = right_maps[equation_index](
                ambient_projections[right_component_indices[equation_index]](
                    ambient_generator,
                ),
            )

            compatibility_value += compatibility_injections[equation_index](
                left_value - right_value,
            )
        end

        push!(images, compatibility_value)
    end

    return _homomorphism_from_generator_images(
        ambient_module,
        compatibility_module,
        images,
    )
end

function MackeyFunctorHomModule(
    domain_mf::MackeyFunctor,
    codomain_mf::MackeyFunctor,
)
    domain_mf.context == codomain_mf.context ||
        throw(ArgumentError("The domain and codomain Mackey functors must have the same Mackey context."))
    coefficient_ring(domain_mf) == coefficient_ring(codomain_mf) ||
        throw(ArgumentError("Mackey functors must be defined over the same coefficient ring."))

    R = coefficient_ring(domain_mf)
    T = elem_type(R)
    component_hom_modules = HomModule{T}[
        HomModule(domain_mf.values[i], codomain_mf.values[i])
        for i in eachindex(domain_mf.context.subgroups)
    ]
    component_modules = AbstractAlgebra.FPModule{T}[
        underlying_module(H)
        for H in component_hom_modules
    ]
    ambient_module, ambient_injections, ambient_projections =
        _direct_sum(component_modules)

    (
        compatibility_hom_modules,
        left_component_indices,
        left_maps,
        right_component_indices,
        right_maps,
    ) = _mackey_functor_hom_compatibility_data(
        domain_mf,
        codomain_mf,
        component_hom_modules,
    )

    compatibility_module, compatibility_injections =
        _compatibility_module_and_injections(R, compatibility_hom_modules)
    compatibility_map = _mackey_functor_hom_compatibility_map(
        ambient_module,
        ambient_projections,
        compatibility_module,
        compatibility_injections,
        left_component_indices,
        left_maps,
        right_component_indices,
        right_maps,
    )

    hom_module, inclusion = kernel(compatibility_map)

    return MackeyFunctorHomModule(
        hom_module,
        domain_mf,
        codomain_mf,
        component_hom_modules,
        ambient_module,
        ambient_injections,
        ambient_projections,
        compatibility_hom_modules,
        compatibility_module,
        compatibility_map,
        inclusion,
    )
end

function HomModule(domain_mf::MackeyFunctor, codomain_mf::MackeyFunctor)
    return MackeyFunctorHomModule(domain_mf, codomain_mf)
end

"""
    underlying_module(H::MackeyFunctorHomModule)

Return the finitely presented module representing Mackey functor
homomorphisms.
"""
underlying_module(H::MackeyFunctorHomModule) = H.H

AbstractAlgebra.base_ring(H::MackeyFunctorHomModule) = base_ring(underlying_module(H))
AbstractAlgebra.number_of_generators(H::MackeyFunctorHomModule) = ngens(underlying_module(H))
AbstractAlgebra.gens(H::MackeyFunctorHomModule) = gens(underlying_module(H))
AbstractAlgebra.gen(H::MackeyFunctorHomModule, i::Int) = gen(underlying_module(H), i)
AbstractAlgebra.relations(H::MackeyFunctorHomModule) = relations(underlying_module(H))

"""
    as_homomorphism(H::MackeyFunctorHomModule, x::FPModuleElem)

Convert an element of the underlying Hom module to the represented Mackey
functor homomorphism.
"""
function as_homomorphism(
    H::MackeyFunctorHomModule{T},
    x::AbstractAlgebra.FPModuleElem{T},
) where T <: RingElement
    parent(x) === underlying_module(H) ||
        throw(ArgumentError("Element does not belong to this Mackey functor Hom module."))

    ambient_element = H.inclusion(x)
    components = Generic.ModuleHomomorphism[
        as_homomorphism(
            H.component_hom_modules[i],
            H.ambient_projections[i](ambient_element),
        )
        for i in eachindex(H.component_hom_modules)
    ]

    return MackeyFunctorHomomorphism(H.domain_mf, H.codomain_mf, components)
end

"""
    as_hom_module_element(H::MackeyFunctorHomModule, f::MackeyFunctorHomomorphism)

Convert a Mackey functor homomorphism `f` into the corresponding element of
the underlying Hom module.
"""
function as_hom_module_element(
    H::MackeyFunctorHomModule{T},
    f::MackeyFunctorHomomorphism,
) where T <: RingElement
    f.domain === H.domain_mf ||
        throw(ArgumentError("The Mackey functor homomorphism has the wrong domain."))
    f.codomain === H.codomain_mf ||
        throw(ArgumentError("The Mackey functor homomorphism has the wrong codomain."))

    ambient_element = zero(H.ambient_module)
    for i in eachindex(H.component_hom_modules)
        component_element = as_hom_module_element(
            H.component_hom_modules[i],
            f.components[i],
        )
        ambient_element += H.ambient_injections[i](component_element)
    end

    has_lift, hom_element = has_preimage_with_preimage(H.inclusion, ambient_element)
    has_lift ||
        throw(ArgumentError("The Mackey functor homomorphism does not commute with the Mackey structure maps."))

    return hom_element
end

function Base.show(io::IO, H::MackeyFunctorHomModule)
    print(io, "Mackey functor Hom module from ")
    print(io, H.domain_mf)
    print(io, " to ")
    print(io, H.codomain_mf)
end
