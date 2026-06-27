function _nonzero_generator_choice(
    mf::MackeyFunctor,
    level_order,
)
    levels =
        level_order === nothing ? eachindex(mf.context.subgroups) :
        level_order isa Function ? level_order(mf) :
        level_order

    for H_index in levels
        checkbounds(mf.context.subgroups, H_index)
        for x in gens(mf.values[H_index])
            iszero(x) || return H_index, x
        end
    end

    return nothing
end

function _is_zero_mackey_functor(mf::MackeyFunctor; level_order=nothing)::Bool
    return _nonzero_generator_choice(mf, level_order) === nothing
end

function _direct_sum_map_to_common_codomain(
    f::MackeyFunctorHomomorphism,
    g::MackeyFunctorHomomorphism,
)::MackeyFunctorHomomorphism
    f.codomain === g.codomain ||
        throw(ArgumentError("Mackey functor homomorphisms must have the same codomain."))

    domain_sum, = direct_sum(f.domain, g.domain)
    components = Generic.ModuleHomomorphism[]
    for H_index in eachindex(f.context.subgroups)
        source = domain_sum.values[H_index]
        target = f.codomain.values[H_index]

        # The direct-sum presentation uses the generators of the left summand
        # followed by the generators of the right summand.  We can therefore
        # build the map out of the direct sum by assigning the old generator
        # images first and the new universal-map generator images second.
        images = elem_type(target)[
            f.components[H_index](x)
            for x in gens(f.domain.values[H_index])
        ]
        append!(
            images,
            elem_type(target)[
                g.components[H_index](x)
                for x in gens(g.domain.values[H_index])
            ],
        )

        push!(
            components,
            _module_homomorphism_from_images(source, target, images),
        )
    end

    return MackeyFunctorHomomorphism(domain_sum, f.codomain, components)
end

"""
    epimorphism_from_free(M::MackeyFunctor; level_order=nothing, verify::Bool=true)

Construct a levelwise-surjective homomorphism from a finite direct sum of free
Mackey functors onto `M`.

The algorithm repeatedly chooses a nonzero generator in the current cokernel,
lifts it back to `M`, adds the universal map that hits that lift, and then
recomputes the cokernel.  The optional `level_order` controls which subgroup
levels are inspected first; pass either a vector of subgroup indices or a
function `C -> indices` that can depend on the current cokernel.
"""
function epimorphism_from_free(
    mf::MackeyFunctor;
    level_order=nothing,
    verify::Bool=true,
)::MackeyFunctorHomomorphism
    if _is_zero_mackey_functor(mf; level_order=level_order)
        return id_homomorphism(mf)
    end

    current_map = nothing
    current_cokernel = mf
    current_projection = id_homomorphism(mf)

    while true
        choice = _nonzero_generator_choice(current_cokernel, level_order)
        choice === nothing && return current_map

        H_index, cokernel_generator = choice
        lift_to_mf = preimage(
            current_projection.components[H_index],
            cokernel_generator,
        )
        new_universal_map = universal_map(
            mf,
            H_index,
            lift_to_mf;
            verify=verify,
        )

        current_map =
            current_map === nothing ? new_universal_map :
            _direct_sum_map_to_common_codomain(current_map, new_universal_map)

        current_cokernel, current_projection = cokernel(current_map)
    end
end
