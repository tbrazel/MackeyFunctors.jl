"""
    constant_mackey_functor(ctx::MackeyContext, M::AbstractAlgebra.FPModule) -> MackeyFunctor


Given a [`MackeyContext`](@ref) and an ``R``-module ``M``, this method outputs the fixed-point Mackey functor for ``M`` with trivial ``G``-action. This is also called the *constant Mackey functor* valued at ``M``.
"""
function constant_mackey_functor(context::MackeyContext, M::AbstractAlgebra.FPModule)
    R = base_ring(M)
    id_hom = identity_homomorphism(M)
    id_iso = identity_isomorphism(M)

    # Every value is M
    values = fill(M, length(context.subgroups))

    # Every restriction is id_M
    restrictions = fill(id_hom, length(context.covers))

    # A transfer M(H)->M(K) is multiplication by [K:H]
    transfers = eltype(restrictions)[
        R(subgroup_inclusion_index(context, cover_index)) * id_hom
        for cover_index in context.covers
    ]

    # Every conjugation is the identity
    conjugations = fill(id_iso, length(context.generators), length(context.subgroups))
    MackeyFunctor(context, values, restrictions, transfers, conjugations)
end
"""
    constant_mackey_functor(ctx::MackeyContext, R::Ring) -> MackeyFunctor

This method can also be fed a context and a ring ``R``, and it will output the constant Mackey functor valued at ``R`` considered as a free rank one ``R``-module.
"""
function constant_mackey_functor(context::MackeyContext, R::Ring)
    M = free_module(R, 1)
    return constant_mackey_functor(context, M)
end

function burnside_transfer(R, M1, M2, cc1, cc2, H2)
    m = zero_matrix(R, rank(M1), rank(M2))
    for (i, II) in enumerate(cc1)
        I = GAP.Globals.Representative(II)
        j = findfirst(==(GAP.Globals.ConjugacyClassSubgroups(H2, I)), cc2)::Int
        m[i, j] = 1
    end
    ModuleHomomorphism(M1, M2, m)
end

function burnside_restriction(R, M1, M2, cc1, cc2, H1, H2)
    m = zero_matrix(R, rank(M2), rank(M1))
    for (j, JJ) in enumerate(cc2)
        J = GAP.Globals.Representative(JJ)
        for (h, n) in GAP.Globals.DoubleCosetRepsAndSizes(H2, H1, J)
            L = GAP.Globals.Intersection(H1, J ^ inv(h))
            i = findfirst(==(GAP.Globals.ConjugacyClassSubgroups(H1, L)), cc1)::Int
            m[j, i] += n * GAP.Globals.Size(L) ÷ (GAP.Globals.Size(H1) * GAP.Globals.Size(J))
        end
    end
    ModuleHomomorphism(M2, M1, m)
end

function burnside_conjugation(R, mc, conj_classes, gi, Hi, M)
    g = mc.generators[gi]
    m = zero_matrix(R, rank(M), rank(M))
    for (j, JJ) in enumerate(conj_classes[Hi])
        L = GAP.Globals.Representative(JJ) ^ inv(g)
        gHi = mc.generator_left_conjugation_matrix[gi, Hi]
        i = findfirst(==(GAP.Globals.ConjugacyClassSubgroups(mc.subgroups[gHi], L)), conj_classes[gHi])::Int
        m[j, i] = 1
    end
    ModuleIsomorphism(M, M, m)
end

"""
    burnside_mackey_functor(mc::MackeyContext, R::Ring = ZZ) -> MackeyFunctor

Return the Burnside Mackey functor for the group specified by `mc` and the coefficient ring `R`.
"""
function burnside_mackey_functor(mc::MackeyContext, R::Ring=ZZ)
    conj_classes = [collect(GAP.Globals.ConjugacyClassesSubgroups(H)) for H in mc.subgroups]
    values = [free_module(R, length(cc)) for cc in conj_classes]

    cover_transfers = Generic.ModuleHomomorphism[
        burnside_transfer(R, values[i], values[j], conj_classes[i], conj_classes[j], mc.subgroups[j])
        for (i, j) in mc.covers
    ]
    cover_restrictions = Generic.ModuleHomomorphism[
        burnside_restriction(R, values[i], values[j], conj_classes[i], conj_classes[j], mc.subgroups[i], mc.subgroups[j])
        for (i, j) in mc.covers
    ]
    generator_conjugations = Generic.ModuleIsomorphism[
        burnside_conjugation(R, mc, conj_classes, gi, Hi, values[Hi])
        for gi in eachindex(mc.generators), Hi in eachindex(mc.subgroups)
    ]

    MackeyFunctor(mc, values, cover_restrictions, cover_transfers, generator_conjugations)  # TODO: turn off argument checking
end

"""
    burnside_mackey_functor(G, R::Ring = ZZ) -> MackeyFunctor

Return the Burnside Mackey functor for the group `G` and the coefficient ring `R`.
"""
burnside_mackey_functor(G::GapObj, R::Ring=ZZ) = burnside_mackey_functor(MackeyContext(G), R)

"""
    free_mackey_functor(mf::MackeyFunctor, i::SubgroupIndex) -> MackeyFunctor

Given a Mackey functor `M` and a subgroup index `i` corresponding to a subgroup ``H\\le G``, this returns the shifted Mackey functor ``M_H``.
"""
function free_mackey_functor(mf::MackeyFunctor, i::SubgroupIndex; verify::Bool=true)
    shift(mf, i, verify)
end

# fixedpoints Mackey functor

function fixedpoint_with_inclusion(gm::GModule, H)
    idM = identity_homomorphism(gm.M)
    subM = sub(gm.M, generators(gm.M))  # TODO: more efficient way?
    foldl(GAP.Globals.GeneratorsOfGroup(H); init=subM) do (N, f), g
        word = generator_word(gm.context, g)
        K, _ = kernel(map_extension(idM, gm.generator_actions, word) - idM)
        N1, f1 = intersect(N, K)
        N1, f1 * f
    end
end

function fixedpoint_transfer(gm, i, j, sub_incl)
    idM = identity_homomorphism(gm.M)
    f = sum(GAP.Globals.RightCosets(gm.context.subgroups[j], gm.context.subgroups[i])) do rc
        h = inv(GAP.Globals.Representative(rc))
        map_extension(idM, gm.generator_actions, generator_word(gm.context, h))
    end
    ModuleHomomorphism(first(sub_incl[i]), first(sub_incl[j]), submodules_matrix(f, sub_incl[i], sub_incl[j]))
end

function fixedpoint_conjugation(gm, gi, Hi, sub_incl)
    j = gm.context.generator_left_conjugation_matrix[gi, Hi]
    ModuleIsomorphism(first(sub_incl[Hi]), first(sub_incl[j]), submodules_matrix(gm.generator_actions[gi], sub_incl[Hi], sub_incl[j]))
end

"""
    fixedpoint_mackey_functor(gm::GModule) -> MackeyFunctor

Return the Mackey functor whose levels are the invariant submodules of `gm` under the subgroups of `gm.context.group`.
"""
function fixedpoint_mackey_functor(gm::GModule)
    sub_incl = [fixedpoint_with_inclusion(gm, H) for H in gm.context.subgroups]

    values = map(first, sub_incl)
    restrictions = [ModuleHomomorphism(values[j], values[i], submodules_matrix(sub_incl[j], sub_incl[i])) for (i, j) in gm.context.covers]
    transfers = [fixedpoint_transfer(gm, i, j, sub_incl) for (i, j) in gm.context.covers]
    conjugations = [fixedpoint_conjugation(gm, gi, Hi, sub_incl)
                    for gi in eachindex(gm.context.generators), Hi in eachindex(gm.context.subgroups)]

    MackeyFunctor(gm.context, values, restrictions, transfers, conjugations)
end


function zero_mackey_functor(context::MackeyContext,R::Ring=ZZ)
    zeromod = free_module(R,0)
    zeromap = ModuleIsomorphism(zeromod, zeromod, ZZ[;])

    values = fill(zeromod, length(context.subgroups))
    cover_restrictions = [as_homomorphism(zeromap) for _ in context.covers]
    cover_transfers = [as_homomorphism(zeromap) for _ in context.covers]
    generator_conjugations = [zeromap for _ in context.generators, _ in context.subgroups]

    MackeyFunctor(context, values, cover_restrictions, cover_transfers, generator_conjugations)
end
