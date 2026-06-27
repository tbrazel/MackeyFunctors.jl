struct ShiftOrbitDecomposition
    # The shifted Mackey functor stores the actual value module.  The
    # decomposition only needs the structure maps identifying its orbit
    # summands inside that value.
    injections::Vector{Generic.ModuleHomomorphism}
    projections::Vector{Generic.ModuleHomomorphism}
end

struct ShiftedMackeyFunctor <: AbstractShiftedMackeyFunctor
    underlying_mackey_functor::MackeyFunctor
    decompositions::Vector{ShiftOrbitDecomposition}
end


function _shift(
    mf::MackeyFunctor,
    H_index::SubgroupIndex;
    verify::Bool=true,
)::ShiftedMackeyFunctor
    # Return the ShiftedMackeyFunctor if it is cached
    if haskey(mf.shift_cache, H_index)
        return mf.shift_cache[H_index]
    end
    ctx = mf.context
    !verify || checkbounds(ctx.subgroups, H_index)

    values = AbstractAlgebra.FPModule[]
    decompositions = ShiftOrbitDecomposition[]
    for K_index in eachindex(ctx.subgroups)
        shifted_value, decomposition =
            _shift_orbit_decomposition(mf, H_index, K_index)
        push!(values, shifted_value)
        push!(decompositions, decomposition)
    end

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
            values,
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
            values,
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

    mf.shift_cache[H_index] = ShiftedMackeyFunctor(MackeyFunctor(
            ctx,
            values,
            cover_restrictions,
            cover_transfers,
            generator_conjugations;
            verify=verify,
        ), decompositions)

    return mf.shift_cache[H_index]
end



function _shift_orbit_decomposition(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    K_index::SubgroupIndex,
)
    ctx = mf.context

    stabilizer_indices = [
        info.left_intersection_conjugated_right_index
        for info in _shift_decomposition_double_coset_infos(ctx, H_index, K_index)
    ]

    summands = AbstractAlgebra.FPModule[
        value(mf, stabilizer_index)
        for stabilizer_index in stabilizer_indices
    ]
    shifted_value, injections, projections = direct_sum(summands)

    return shifted_value, ShiftOrbitDecomposition(
        Generic.ModuleHomomorphism[injections...],
        Generic.ModuleHomomorphism[projections...],
    )
end

function _shift_product_maps(
    mf::MackeyFunctor,
    H_index::SubgroupIndex,
    values::Vector{AbstractAlgebra.FPModule},
    decompositions::Vector{ShiftOrbitDecomposition},
    source_index::SubgroupIndex,
    target_index::SubgroupIndex,
    transporter::GroupElement,
    verify::Bool,
)
    ctx = mf.context
    H = ctx.subgroups[H_index]
    target_subgroup = ctx.subgroups[target_index]

    double_coset_source =
        _shift_decomposition_double_coset_infos(ctx, H_index, source_index)
    double_coset_target =
        _shift_decomposition_double_coset_infos(ctx, H_index, target_index)
    source_decomposition = decompositions[source_index]
    target_decomposition = decompositions[target_index]

    # Start with zero maps between the two direct sums, then add one summand
    # map at a time using the injections and projections returned by
    # direct_sum.
    # If f_ij : source_summand_i -> target_summand_j is the map attached to one
    # product-orbit component, then
    #
    #     source_projection_i * f_ij * target_injection_j
    #
    # is the corresponding map from the whole source direct sum to the whole
    # target direct sum.
    covariant_map = zero_homomorphism(
        values[source_index],
        values[target_index],
    )
    contravariant_map = zero_homomorphism(
        values[target_index],
        values[source_index],
    )

    for source_orbit_index in eachindex(double_coset_source)

        # x is a double coset representative for H\G/K_source
        x = double_coset_source[source_orbit_index].representative

        # source_stabilizer_index is the index of H\cap x K_source x^-1
        source_stabilizer_index =
            double_coset_source[
                source_orbit_index
            ].left_intersection_conjugated_right_index

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
            double_coset_target,
            image_representative,
        )

        target_stabilizer_index =
            double_coset_target[
                target_orbit_index
            ].left_intersection_conjugated_right_index
        target_stabilizer = ctx.subgroups[target_stabilizer_index]

        covariant_block, contravariant_block = _shift_orbit_map_blocks(
            mf,
            source_stabilizer_index,
            target_stabilizer_index,
            target_stabilizer,
            left_transporter,
            verify,
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

function _shift_orbit_map_blocks(
    mf::MackeyFunctor,
    source_stabilizer_index::SubgroupIndex,
    target_stabilizer_index::SubgroupIndex,
    target_stabilizer::Group,
    orbit_transporter::GroupElement,
    verify::Bool,
)
    ctx = mf.context

    # This helper is the small piece of Mackey-functor algebra shared by all
    # product maps used in shifts.
    #
    # Suppose one transitive source orbit has stabilizer S and the target orbit
    # chosen in our double-coset decomposition has stabilizer T.  A G-map between
    # these two transitive orbits is determined by an element b of G:
    #
    #     alpha : G/S -> G/T,     qS |-> q b T.
    #
    # The formula is well-defined exactly when S <= bTb^-1.  Write
    #
    #     U = bTb^-1.
    #
    # Then alpha factors as
    #
    #     G/S --projection--> G/U --isomorphism--> G/T,
    #     qS |-> qU              qU |-> q b T.
    #
    # A Mackey functor is covariant for projections and contravariant for the
    # same projections, so this one orbit map gives two module maps:
    #
    #   * covariant:     M(S) -> M(U) -> M(T),
    #                    transfer along S <= U, then conjugation by b^-1;
    #
    #   * contravariant: M(T) -> M(U) -> M(S),
    #                    conjugation by b, then restriction along S <= U.
    #
    # GAP writes K^g for g^-1*K*g, so bTb^-1 is `T^(b^-1)`.
    conjugated_target_stabilizer =
        target_stabilizer^(orbit_transporter^-1)

    conjugated_target_stabilizer_index =
        subgroup_index(ctx, conjugated_target_stabilizer)

    if verify
        is_subgroup(
            ctx,
            source_stabilizer_index,
            conjugated_target_stabilizer_index,
        ) || throw(ArgumentError("Product orbit map did not preserve stabilizers."))
    end

    covariant_block =
        transfer(
            mf,
            source_stabilizer_index,
            conjugated_target_stabilizer_index,
        ) * conjugation(
            mf,
            orbit_transporter^-1,
            conjugated_target_stabilizer_index,
        )

    contravariant_block =
        conjugation(
            mf,
            orbit_transporter,
            target_stabilizer_index,
        ) * restriction(
            mf,
            source_stabilizer_index,
            conjugated_target_stabilizer_index,
        )

    return covariant_block, contravariant_block
end

function _shift_first_factor_product_maps(
    mf::MackeyFunctor,
    source_shift::ShiftedMackeyFunctor,
    target_shift::ShiftedMackeyFunctor,
    source_H_index::SubgroupIndex,
    target_H_index::SubgroupIndex,
    K_index::SubgroupIndex,
    transporter::GroupElement,
    verify::Bool,
)
    ctx = mf.context
    target_H = ctx.subgroups[target_H_index]
    K = ctx.subgroups[K_index]

    double_coset_source =
        _shift_decomposition_double_coset_infos(ctx, source_H_index, K_index)
    double_coset_target =
        _shift_decomposition_double_coset_infos(ctx, target_H_index, K_index)
    source_decomposition = source_shift.decompositions[K_index]
    target_decomposition = target_shift.decompositions[K_index]

    source_value = source_shift.underlying_mackey_functor.values[K_index]
    target_value = target_shift.underlying_mackey_functor.values[K_index]

    covariant_map = zero_homomorphism(source_value, target_value)
    contravariant_map = zero_homomorphism(target_value, source_value)

    for source_orbit_index in eachindex(double_coset_source)
        # The source summand indexed by x represents the transitive orbit of
        #
        #     (H_s, xK) in (G/H_s) x (G/K).
        #
        # Its stabilizer is S = H_s ∩ xKx^-1.
        x = double_coset_source[source_orbit_index].representative
        source_stabilizer_index =
            double_coset_source[
                source_orbit_index
            ].left_intersection_conjugated_right_index

        # We are handling first-coordinate orbit maps
        #
        #     G/H_s -> G/H_t,     qH_s |-> q r H_t.
        #
        # At level K, the source representative (H_s, xK) maps to
        #
        #     (rH_t, xK).
        #
        # Our target decomposition, however, uses representatives of the form
        # (H_t, yK).  Acting by r^-1 moves the first coordinate back to H_t, so
        # the target double coset is determined by r^-1*x in H_t\G/K.
        image_representative = transporter^-1 * x
        target_orbit_index, left_transporter = _find_shift_target_orbit(
            target_H,
            K,
            double_coset_target,
            image_representative,
        )

        target_stabilizer_index =
            double_coset_target[
                target_orbit_index
            ].left_intersection_conjugated_right_index
        target_stabilizer = ctx.subgroups[target_stabilizer_index]

        # If r^-1*xK = a*yK with a in H_t, then
        #
        #     xK = r*a*yK.
        #
        # Thus the image point (rH_t, xK) is obtained from the chosen target
        # representative (H_t, yK) by the single group element b = r*a.  The
        # induced map on this pair of transitive orbits is therefore
        #
        #     G/S -> G/T,     qS |-> q b T.
        orbit_transporter = transporter * left_transporter
        covariant_block, contravariant_block = _shift_orbit_map_blocks(
            mf,
            source_stabilizer_index,
            target_stabilizer_index,
            target_stabilizer,
            orbit_transporter,
            verify,
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
    target_infos::Vector{DoubleCosetInfo},
    image_representative::GroupElement,
)::Tuple{Int,GroupElement}
    for (target_orbit_index, info) in enumerate(target_infos)
        y = info.representative
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

function _check_shift_orbit_map(
    ctx::MackeyContext,
    source_H_index::SubgroupIndex,
    target_H_index::SubgroupIndex,
    transporter::GroupElement,
)
    checkbounds(ctx.subgroups, source_H_index)
    checkbounds(ctx.subgroups, target_H_index)

    # A formula qH_s |-> q r H_t gives a well-defined G-map G/H_s -> G/H_t
    # exactly when H_s is contained in r H_t r^-1.
    source_H = ctx.subgroups[source_H_index]
    target_H = ctx.subgroups[target_H_index]
    conjugated_target_H = target_H^(transporter^-1)

    Bool(GAP.Globals.IsSubgroup(conjugated_target_H, source_H)) ||
        throw(ArgumentError("The orbit formula qH_source |-> q*r*H_target is not well-defined."))

    return nothing
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
function shift(phi::MackeyFunctorHomomorphism, H_index::SubgroupIndex)::MackeyFunctorHomomorphism
    ctx = phi.context

    # Domain/codomain as ShiftedMackeyFunctor types - will be pulled from cache if they were already built
    new_domain = _shift(phi.domain, H_index)
    new_codomain = _shift(phi.codomain, H_index)

    values = Vector{Generic.ModuleHomomorphism}()

    # Build the values of phi_H
    for k in eachindex(ctx.subgroups)
        # Start with the zero map
        value = zero_homomorphism(new_domain.underlying_mackey_functor.values[k], new_codomain.underlying_mackey_functor.values[k])

        # We write G = \cup_x HxK.  The component indexed by x is the
        # stabilizer H ∩ xKx^-1 stored in the shared double-coset info.
        for (i, info) in enumerate(
            _shift_decomposition_double_coset_infos(ctx, H_index, k),
        )
            stabilizer_index = info.left_intersection_conjugated_right_index
            value +=
                new_domain.decompositions[k].projections[i] *
                phi.components[stabilizer_index] *
                new_codomain.decompositions[k].injections[i]
        end
        push!(values, value)
    end

    return MackeyFunctorHomomorphism(
        new_domain.underlying_mackey_functor,
        new_codomain.underlying_mackey_functor,
        values
    )
end

"""
    shift_transfer(M::MackeyFunctor, H1_index, H2_index; verify::Bool=true)

Return the canonical homomorphism ``M_{H_1} \\to M_{H_2}`` induced by the
projection ``G/H_1 \\to G/H_2`` when ``H_1 \\le H_2``.
"""
function shift_transfer(
    mf::MackeyFunctor,
    H1_index::SubgroupIndex,
    H2_index::SubgroupIndex;
    verify::Bool=true,
)::MackeyFunctorHomomorphism
    ctx = mf.context
    checkbounds(ctx.subgroups, H1_index)
    checkbounds(ctx.subgroups, H2_index)
    is_subgroup(ctx, H1_index, H2_index) ||
        throw(ArgumentError("The first subgroup must be contained in the second subgroup."))

    source_shift = _shift(mf, H1_index; verify=verify)
    target_shift = _shift(mf, H2_index; verify=verify)
    identity_element = GAP.Globals.One(ctx.group)

    components = Generic.ModuleHomomorphism[]
    for K_index in eachindex(ctx.subgroups)
        covariant_map, = _shift_first_factor_product_maps(
            mf,
            source_shift,
            target_shift,
            H1_index,
            H2_index,
            K_index,
            identity_element,
            verify,
        )
        push!(components, covariant_map)
    end

    return MackeyFunctorHomomorphism(
        source_shift.underlying_mackey_functor,
        target_shift.underlying_mackey_functor,
        components,
    )
end

"""
    shift_restriction(M::MackeyFunctor, H1_index, H2_index; verify::Bool=true)

Return the canonical homomorphism ``M_{H_2} \\to M_{H_1}`` induced
contravariantly by the projection ``G/H_1 \\to G/H_2`` when ``H_1 \\le H_2``.
"""
function shift_restriction(
    mf::MackeyFunctor,
    H1_index::SubgroupIndex,
    H2_index::SubgroupIndex;
    verify::Bool=true,
)::MackeyFunctorHomomorphism
    ctx = mf.context
    checkbounds(ctx.subgroups, H1_index)
    checkbounds(ctx.subgroups, H2_index)
    is_subgroup(ctx, H1_index, H2_index) ||
        throw(ArgumentError("The first subgroup must be contained in the second subgroup."))

    source_shift = _shift(mf, H1_index; verify=verify)
    target_shift = _shift(mf, H2_index; verify=verify)
    identity_element = GAP.Globals.One(ctx.group)

    components = Generic.ModuleHomomorphism[]
    for K_index in eachindex(ctx.subgroups)
        _, contravariant_map = _shift_first_factor_product_maps(
            mf,
            source_shift,
            target_shift,
            H1_index,
            H2_index,
            K_index,
            identity_element,
            verify,
        )
        push!(components, contravariant_map)
    end

    return MackeyFunctorHomomorphism(
        target_shift.underlying_mackey_functor,
        source_shift.underlying_mackey_functor,
        components,
    )
end

"""
    shift_conjugation(M::MackeyFunctor, g::GroupElement, H_index; verify::Bool=true)

Return the canonical homomorphism ``M_H \\to M_{gHg^{-1}}`` induced by the
orbit isomorphism ``G/H \\to G/(gHg^{-1})`` given by
``qH \\mapsto qg^{-1}(gHg^{-1})``.
"""
function shift_conjugation(
    mf::MackeyFunctor,
    g::GroupElement,
    H_index::SubgroupIndex;
    verify::Bool=true,
)::MackeyFunctorHomomorphism
    ctx = mf.context
    checkbounds(ctx.subgroups, H_index)

    target_H = ctx.subgroups[H_index]^(g^-1)
    target_H_index = subgroup_index(ctx, target_H)
    source_shift = _shift(mf, H_index; verify=verify)
    target_shift = _shift(mf, target_H_index; verify=verify)
    transporter = g^-1

    if verify
        _check_shift_orbit_map(ctx, H_index, target_H_index, transporter)
    end

    components = Generic.ModuleHomomorphism[]
    for K_index in eachindex(ctx.subgroups)
        covariant_map, = _shift_first_factor_product_maps(
            mf,
            source_shift,
            target_shift,
            H_index,
            target_H_index,
            K_index,
            transporter,
            verify,
        )
        push!(components, covariant_map)
    end

    return MackeyFunctorHomomorphism(
        source_shift.underlying_mackey_functor,
        target_shift.underlying_mackey_functor,
        components,
    )
end
