"""
    free_mackey_functor(ctx::MackeyContext, i::SubgroupIndex, R::Ring = ZZ) -> MackeyFunctor

Return the free Mackey functor at the subgroup indexed by `i`.
"""
function free_mackey_functor(
    context::MackeyContext,
    i::SubgroupIndex,
    R::Ring=ZZ;
    verify::Bool=true,
)
    return shift(burnside_mackey_functor(context, R), i; verify=verify)
end

function _shifted_burnside_mackey_functor(
    context::MackeyContext,
    H_index::SubgroupIndex,
    R::Ring;
    verify::Bool=true,
)
    checkbounds(context.subgroups, H_index)

    conj_classes = _burnside_conjugacy_classes(context)
    burnside = _burnside_mackey_functor(context, R, conj_classes)
    shifted_burnside = _shift(burnside, H_index; verify=verify)
    basis_subgroup_indices =
        _burnside_basis_subgroup_indices(context, conj_classes)

    return burnside, shifted_burnside, basis_subgroup_indices
end

function _burnside_basis_subgroup_indices(
    context::MackeyContext,
    conj_classes,
)::Vector{Vector{SubgroupIndex}}
    return [
        SubgroupIndex[
            subgroup_index(context, GAP.Globals.Representative(conjugacy_class))
            for conjugacy_class in subgroup_conj_classes
        ]
        for subgroup_conj_classes in conj_classes
    ]
end

function _burnside_basis_subgroup_indices(
    context::MackeyContext,
)::Vector{Vector{SubgroupIndex}}
    return _burnside_basis_subgroup_indices(
        context,
        _burnside_conjugacy_classes(context),
    )
end

function _identity_double_coset_orbit_index(
    context::MackeyContext,
    H_index::SubgroupIndex,
)::Int
    H = context.subgroups[H_index]
    identity_element = GAP.Globals.One(context.group)
    orbit_index = findfirst(
        info -> Bool(GAP.Globals.IN(
            identity_element,
            GAP.Globals.DoubleCoset(H, info.representative, H),
        )),
        _shift_decomposition_double_coset_infos(context, H_index, H_index),
    )

    orbit_index === nothing &&
        throw(ArgumentError("Could not locate the identity orbit in G/H x G/H."))

    return orbit_index
end

function _universal_element(
    burnside::MackeyFunctor,
    shifted_burnside::ShiftedMackeyFunctor,
    H_index::SubgroupIndex,
    basis_subgroup_indices::Vector{Vector{SubgroupIndex}},
)
    orbit_index = _identity_double_coset_orbit_index(burnside.context, H_index)
    basis_index = findfirst(==(H_index), basis_subgroup_indices[H_index])
    basis_index === nothing &&
        throw(ArgumentError("Could not locate the H/H basis element in A(H)."))

    return shifted_burnside.decompositions[H_index].injections[orbit_index](
        gen(burnside.values[H_index], basis_index),
    )
end

"""
    universal_element(ctx::MackeyContext, H_index::SubgroupIndex, R::Ring=ZZ; verify::Bool=true)

Return `(A_H, u_H)`, where `A_H` is the free Mackey functor at the
subgroup indexed by `H_index` and `u_H in A_H(H)` is the element represented
by the identity span `G/H <- G/H -> G/H`.
"""
function universal_element(
    context::MackeyContext,
    H_index::SubgroupIndex,
    R::Ring=ZZ;
    verify::Bool=true,
)
    burnside, shifted_burnside, basis_subgroup_indices =
        _shifted_burnside_mackey_functor(context, H_index, R; verify=verify)

    return (
        shifted_burnside.underlying_mackey_functor,
        _universal_element(
            burnside,
            shifted_burnside,
            H_index,
            basis_subgroup_indices,
        ),
    )
end

function _universal_span_action(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    L_index::SubgroupIndex,
    conjugated_K_index::SubgroupIndex,
    representative::GroupElement,
    x::AbstractAlgebra.FPModuleElem,
)
    # A basis element in the summand for a product orbit is a span
    #
    #     G/H <- G/L -> G/K,
    #
    # where L <= H cap rKr^-1 and the right map is qL |-> q r K.
    # Acting on x in M(H) therefore means restrict to L, transfer to rKr^-1,
    # then conjugate by r^-1 to land at K.
    return (
        restriction(mf, L_index, H_index) *
        transfer(mf, L_index, conjugated_K_index) *
        conjugation(mf, representative^-1, conjugated_K_index)
    )(x)
end

function _universal_map_component(
    mf::MackeyFunctor,
    burnside::MackeyFunctor,
    shifted_burnside::ShiftedMackeyFunctor,
    H_index::SubgroupIndex,
    K_index::SubgroupIndex,
    x::AbstractAlgebra.FPModuleElem,
    basis_subgroup_indices::Vector{Vector{SubgroupIndex}},
)
    context = mf.context
    source = shifted_burnside.underlying_mackey_functor.values[K_index]
    target = mf.values[K_index]
    component = zero_homomorphism(source, target)
    decomposition = shifted_burnside.decompositions[K_index]

    for (orbit_index, info) in enumerate(
        _shift_decomposition_double_coset_infos(context, H_index, K_index),
    )
        stabilizer_index = info.left_intersection_conjugated_right_index
        conjugated_K_index = subgroup_index(
            context,
            context.subgroups[K_index]^(info.representative^-1),
        )
        images = elem_type(target)[
            _universal_span_action(
                mf,
                H_index,
                L_index,
                conjugated_K_index,
                info.representative,
                x,
            )
            for L_index in basis_subgroup_indices[stabilizer_index]
        ]
        summand_map = ModuleHomomorphism(
            burnside.values[stabilizer_index],
            target,
            images,
        )

        component += decomposition.projections[orbit_index] * summand_map
    end

    return component
end

function _universal_map(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    x::AbstractAlgebra.FPModuleElem,
    burnside::MackeyFunctor,
    shifted_burnside::ShiftedMackeyFunctor,
    basis_subgroup_indices::Vector{Vector{SubgroupIndex}},
)::MackeyFunctorHomomorphism
    components = Generic.ModuleHomomorphism[
        _universal_map_component(
            mf,
            burnside,
            shifted_burnside,
            H_index,
            K_index,
            x,
            basis_subgroup_indices,
        )
        for K_index in eachindex(mf.context.subgroups)
    ]

    return MackeyFunctorHomomorphism(
        shifted_burnside.underlying_mackey_functor,
        mf,
        components,
    )
end

"""
    universal_map(M::MackeyFunctor, H_index::SubgroupIndex, x; verify::Bool=true)

Return the Mackey functor homomorphism `A_H -> M` that sends the universal
element of the free Mackey functor `A_H` to `x in M(H)`.
"""
function universal_map(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    x::AbstractAlgebra.FPModuleElem;
    verify::Bool=true,
)::MackeyFunctorHomomorphism
    checkbounds(mf.context.subgroups, H_index)
    parent(x) === mf.values[H_index] ||
        throw(ArgumentError("The universal-map element must belong to M(H_index)."))

    burnside, shifted_burnside, basis_subgroup_indices =
        _shifted_burnside_mackey_functor(
            mf.context,
            H_index,
            coefficient_ring(mf);
            verify=verify,
        )

    return _universal_map(
        mf,
        H_index,
        x,
        burnside,
        shifted_burnside,
        basis_subgroup_indices,
    )
end
