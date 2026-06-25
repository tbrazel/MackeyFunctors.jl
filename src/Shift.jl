struct ShiftOrbitDecomposition
    representatives::Vector{GroupElement}
    stabilizer_indices::Vector{SubgroupIndex}
    value::AbstractAlgebra.FPModule
    injections::Vector{Generic.ModuleHomomorphism}
    projections::Vector{Generic.ModuleHomomorphism}
end

function _zero_direct_sum_relation_row(R::Ring, n::Int)
    return [zero(R) for _ in 1:n]
end

function _homomorphism_from_generator_images(
    domain::AbstractAlgebra.FPModule{T},
    codomain::AbstractAlgebra.FPModule{T},
    images::Vector{<:AbstractAlgebra.FPModuleElem{T}},
) where T <: RingElement
    if ngens(domain) == 0
        return ModuleHomomorphism(
            domain,
            codomain,
            zero_matrix(base_ring(domain), 0, ngens(codomain)),
        )
    end

    return ModuleHomomorphism(domain, codomain, images)
end

function _finite_direct_sum(summands::AbstractVector{<:AbstractAlgebra.FPModule})
    isempty(summands) && throw(ArgumentError("Cannot take the direct sum of no modules."))

    R = base_ring(first(summands))
    all(summand -> base_ring(summand) == R, summands) ||
        throw(ArgumentError("Direct-sum summands must have the same base ring."))

    T = elem_type(R)
    return _finite_direct_sum(AbstractAlgebra.FPModule{T}[summands...])
end

function _finite_direct_sum(
    summands::Vector{<:AbstractAlgebra.FPModule{T}},
) where T <: RingElement
    isempty(summands) && throw(ArgumentError("Cannot take the direct sum of no modules."))

    R = base_ring(first(summands))
    all(summand -> base_ring(summand) == R, summands) ||
        throw(ArgumentError("Direct-sum summands must have the same base ring."))

    offsets = Int[]
    total_generators = 0
    for summand in summands
        push!(offsets, total_generators)
        total_generators += ngens(summand)
    end

    free_sum = free_module(R, total_generators)
    relation_generators = elem_type(free_sum)[]
    for (summand_index, summand) in enumerate(summands)
        offset = offsets[summand_index]

        # The direct-sum presentation is formed by placing each summand's
        # relations in its own coordinate block.
        for relation in relations(summand)
            row = _zero_direct_sum_relation_row(R, total_generators)
            for generator_index in 1:ngens(summand)
                row[offset + generator_index] = relation[1, generator_index]
            end
            push!(relation_generators, free_sum(row))
        end
    end

    if isempty(relation_generators)
        direct_sum_module = free_sum
        projection_from_free = identity_homomorphism(free_sum)
    else
        relation_submodule, = sub(free_sum, relation_generators)
        direct_sum_module, projection_from_free = quo(free_sum, relation_submodule)
    end

    injections = Vector{Generic.ModuleHomomorphism{T}}(undef, length(summands))
    projections = Vector{Generic.ModuleHomomorphism{T}}(undef, length(summands))

    for (summand_index, summand) in enumerate(summands)
        offset = offsets[summand_index]

        injection_images = elem_type(direct_sum_module)[
            projection_from_free(gen(free_sum, offset + generator_index))
            for generator_index in 1:ngens(summand)
        ]
        injections[summand_index] = _homomorphism_from_generator_images(
            summand,
            direct_sum_module,
            injection_images,
        )

        projection_images = elem_type(summand)[]
        for direct_sum_generator in gens(direct_sum_module)
            lift = preimage(projection_from_free, direct_sum_generator)
            push!(
                projection_images,
                summand([
                    lift[offset + generator_index]
                    for generator_index in 1:ngens(summand)
                ]),
            )
        end
        projections[summand_index] = _homomorphism_from_generator_images(
            direct_sum_module,
            summand,
            projection_images,
        )
    end

    return direct_sum_module, injections, projections
end

"""
    shift(M::MackeyFunctor, H_index::SubgroupIndex; verify::Bool=true)

Return the ``G/H``-shift of `M`.

At a subgroup ``K`` this Mackey functor has value
```math
M_{G/H}(G/K) = M((G/H) \\times (G/K)).
```
The implementation decomposes the product into transitive ``G``-orbits and
then uses the existing restriction, transfer, and conjugation maps of `M` on
each orbit summand.
"""
function shift(
    mf::MackeyFunctor,
    H_index::SubgroupIndex;
    verify::Bool=true,
)::MackeyFunctor
    ctx = mf.context
    !verify || checkbounds(ctx.subgroups, H_index)

    decompositions = [
        _shift_orbit_decomposition(mf, H_index, K_index)
        for K_index in eachindex(ctx.subgroups)
    ]

    values = AbstractAlgebra.FPModule[
        decomposition.value
        for decomposition in decompositions
    ]

    identity_element = GAP.Globals.One(ctx.group)
    cover_restrictions = Generic.ModuleHomomorphism[]
    cover_transfers = Generic.ModuleHomomorphism[]

    for (K_index, L_index) in ctx.covers
        # The cover K <= L gives a G-map G/K -> G/L by qK |-> qL.
        # On products this is
        #
        #     (G/H) x (G/K) -> (G/H) x (G/L).
        #
        # The shifted transfer is the covariant map along this product map, and
        # the shifted restriction is the contravariant map.
        covariant_map, contravariant_map = _shift_product_maps(
            mf,
            H_index,
            decompositions,
            K_index,
            L_index,
            identity_element,
            verify,
        )

        push!(cover_transfers, covariant_map)
        push!(cover_restrictions, contravariant_map)
    end

    generator_conjugations = Matrix{Generic.ModuleIsomorphism}(
        undef,
        length(ctx.generators),
        length(ctx.subgroups),
    )

    for generator_index in eachindex(ctx.generators), K_index in eachindex(ctx.subgroups)
        g = ctx.generators[generator_index]
        gKginv_index = ctx.generator_left_conjugation_matrix[generator_index, K_index]

        # For an ordinary Mackey functor, the generator conjugation c_g at K is
        # the contravariant map attached to the G-orbit isomorphism
        #
        #     G/(gKg^-1) -> G/K,     q(gKg^-1) |-> qgK.
        #
        # Applying the shift means we multiply this orbit isomorphism by G/H:
        #
        #     (G/H) x G/(gKg^-1) -> (G/H) x G/K,
        #     (pH, q(gKg^-1)) |-> (pH, qgK).
        #
        # Since a Mackey functor is contravariant on restrictions/conjugations,
        # this product isomorphism induces
        #
        #     M_{G/H}(K) -> M_{G/H}(gKg^-1),
        #
        # which is exactly the generator conjugation map needed for the shifted
        # functor.  The helper below constructs both variance directions for a
        # product map whose second-coordinate formula is qK_source |-> q r K_target.
        # Here K_source = gKg^-1, K_target = K, and r = g, so we keep the
        # contravariant direction.
        _, contravariant_map = _shift_product_maps(
            mf,
            H_index,
            decompositions,
            gKginv_index,
            K_index,
            g,
            verify,
        )

        generator_conjugations[generator_index, K_index] = ModuleIsomorphism(
            values[K_index],
            values[gKginv_index],
            matrix(contravariant_map),
        )
    end

    return MackeyFunctor(
        ctx,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations;
        verify=verify,
    )
end

function _shift_orbit_decomposition(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    K_index::SubgroupIndex,
)::ShiftOrbitDecomposition
    ctx = mf.context
    H = ctx.subgroups[H_index]
    K = ctx.subgroups[K_index]

    # A point of (G/H) x (G/K) can be moved into the form (H, xK).
    # Two such points (H, xK) and (H, yK) are in the same orbit exactly
    # when x and y represent the same double coset in H \ G / K.
    representatives = GroupElement[
        entry[1]
        for entry in GAP.Globals.DoubleCosetRepsAndSizes(ctx.group, H, K)
    ]

    # The stabilizer of (H, xK) is
    #
    #     H cap xKx^-1.
    #
    # GAP writes K^a for a^-1 K a, so xKx^-1 is K^(x^-1).
    stabilizer_indices = SubgroupIndex[
        subgroup_index(
            ctx,
            GAP.Globals.Intersection(H, K^(x^-1)),
        )
        for x in representatives
    ]

    summands = AbstractAlgebra.FPModule[
        value(mf, stabilizer_index)
        for stabilizer_index in stabilizer_indices
    ]
    shifted_value, injections, projections = _finite_direct_sum(summands)

    return ShiftOrbitDecomposition(
        representatives,
        stabilizer_indices,
        shifted_value,
        Generic.ModuleHomomorphism[injections...],
        Generic.ModuleHomomorphism[projections...],
    )
end

function _shift_product_maps(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    decompositions::Vector{ShiftOrbitDecomposition},
    source_index::SubgroupIndex,
    target_index::SubgroupIndex,
    transporter::GroupElement,
    verify::Bool,
)
    ctx = mf.context
    H = ctx.subgroups[H_index]
    target_subgroup = ctx.subgroups[target_index]

    source_decomposition = decompositions[source_index]
    target_decomposition = decompositions[target_index]

    # Start with zero maps between the two direct sums, then add one summand
    # map at a time using the injections and projections returned by
    # _finite_direct_sum.
    # If f_ij : source_summand_i -> target_summand_j is the map attached to one
    # product-orbit component, then
    #
    #     source_projection_i * f_ij * target_injection_j
    #
    # is the corresponding map from the whole source direct sum to the whole
    # target direct sum.
    covariant_map = zero_homomorphism(
        source_decomposition.value,
        target_decomposition.value,
    )
    contravariant_map = zero_homomorphism(
        target_decomposition.value,
        source_decomposition.value,
    )

    for source_orbit_index in eachindex(source_decomposition.representatives)

        # x is a double coset representative for H\G/K_source
        x = source_decomposition.representatives[source_orbit_index]

        # source_stabilizer_index is the index of H\cap x K_source x^-1
        source_stabilizer_index =
            source_decomposition.stabilizer_indices[source_orbit_index]

        # The general product map handled here is induced by
        #
        #     G/K_source -> G/K_target,     qK_source |-> q r K_target.
        #
        # For restrictions/transfers along K <= L, r is the identity.  For
        # generator conjugation by g, r is g and K_source is gKg^-1.
        #
        # Under this map the orbit representative (H, xK_source) lands at
        # (H, xrK_target).  We locate the unique target double-coset
        # representative y with xr in H y K_target, and an element a in H with
        # xrK_target = a y K_target.
        image_representative = x * transporter
        target_orbit_index, left_transporter = _find_shift_target_orbit(
            H,
            target_subgroup,
            target_decomposition.representatives,
            image_representative,
        )

        target_stabilizer_index =
            target_decomposition.stabilizer_indices[target_orbit_index]
        target_stabilizer = ctx.subgroups[target_stabilizer_index]

        # With y and a as above, the restricted product map on this orbit is
        #
        #     G/S -> G/T,     qS |-> q a T,
        #
        # where S is the source stabilizer and T is the target stabilizer.  This
        # factors as
        #
        #     G/S -> G/(aTa^-1) -> G/T.
        #
        # The first arrow is the projection attached to S <= aTa^-1.  The second
        # arrow is the orbit isomorphism whose contravariant map is the
        # conjugation c_a : M(T) -> M(aTa^-1).
        conjugated_target_stabilizer =
            target_stabilizer^(left_transporter^-1)

        # index of aTa^-1
        conjugated_target_stabilizer_index =
            subgroup_index(ctx, conjugated_target_stabilizer)

        if verify
            is_subgroup(
                ctx,
                source_stabilizer_index,
                conjugated_target_stabilizer_index,
            ) || throw(ArgumentError("Product orbit map did not preserve stabilizers."))
        end

        # Let U = aTa^-1.  The product orbit map restricted to this summand is
        #
        #     alpha : G/S -> G/T,     qS |-> q a T.
        #
        # We use the factorization
        #
        #     G/S --projection--> G/U --isomorphism--> G/T,
        #     qS |-> qU              qU |-> q a T.
        #
        # The covariant map M(S) -> M(T) is therefore transfer along S <= U,
        # followed by the covariant direction of the isomorphism G/U -> G/T.
        # The contravariant direction of that isomorphism is c_a : M(T) -> M(U),
        # so its covariant direction is c_{a^-1} : M(U) -> M(T).
        covariant_block =
            transfer(
                mf,
                source_stabilizer_index,
                conjugated_target_stabilizer_index,
            ) * conjugation(
                mf,
                left_transporter^-1,
                conjugated_target_stabilizer_index,
            )

        # The contravariant map M(T) -> M(S) goes through the same factorization
        # in reverse: first c_a : M(T) -> M(U), then restriction along S <= U.
        contravariant_block =
            conjugation(
                mf,
                left_transporter,
                target_stabilizer_index,
            ) * restriction(
                mf,
                source_stabilizer_index,
                conjugated_target_stabilizer_index,
            )

        covariant_map +=
            source_decomposition.projections[source_orbit_index] *
            covariant_block *
            target_decomposition.injections[target_orbit_index]
        contravariant_map +=
            target_decomposition.projections[target_orbit_index] *
            contravariant_block *
            source_decomposition.injections[source_orbit_index]
    end

    return covariant_map, contravariant_map
end

# Input: H, K, some representatives y_1..y_n for double cosets of H\G/K, and some element x which we think about as representing a double coset HxK
# Output: (i,g) such that Hy_iK = HxK and g in H satisfies g y_i K = x K.
function _find_shift_target_orbit(
    H::Group,
    target_subgroup::Group,
    target_representatives::Vector{GroupElement},
    image_representative::GroupElement,
)::Tuple{Int,GroupElement}
    for (target_orbit_index, y) in enumerate(target_representatives)
        Bool(
            GAP.Globals.IN(
                image_representative,
                GAP.Globals.DoubleCoset(H, y, target_subgroup),
            ),
        ) || continue

        return (
            target_orbit_index,
            _find_left_transporter(
                H,
                target_subgroup,
                y,
                image_representative,
            ),
        )
    end

    throw(ArgumentError("Could not locate target product orbit."))
end

function _find_left_transporter(
    H::Group,
    target_subgroup::Group,
    target_representative::GroupElement,
    image_representative::GroupElement,
)::GroupElement
    # At this point image_representative is known to lie in
    #
    #     H * target_representative * target_subgroup.
    #
    # We need an element a in H such that
    #
    #     image_representative * target_subgroup
    #       = a * target_representative * target_subgroup.
    #
    # Equivalently, a lies in
    #
    #     H \cap image_representative * target_subgroup * target_representative^-1.
    #
    # GAP's RightCoset(K, r) represents K*r, so we rewrite this translated
    # subgroup as a right coset:
    #
    #     image_representative * K * target_representative^-1
    #       = (image_representative * K * image_representative^-1)
    #         * (image_representative * target_representative^-1).
    #
    # Since GAP writes K^g for g^-1*K*g, the conjugated subgroup on the right is
    # target_subgroup^(image_representative^-1).  Asking GAP for a representative
    # of the intersection avoids scanning every element of H in Julia.
    transporter_coset = GAP.Globals.RightCoset(
        target_subgroup^(image_representative^-1),
        image_representative * target_representative^-1,
    )
    return GAP.Globals.Representative(
        GAP.Globals.Intersection(H, transporter_coset),
    )
end
