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
            block_homomorphism([current_map; new_universal_map])

        current_cokernel, current_projection = cokernel(current_map)
    end
end
