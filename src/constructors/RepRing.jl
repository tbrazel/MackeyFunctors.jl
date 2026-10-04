# Given:
#       a list of irred. characters [ψ1,..,ψn] of a subgp H
#       a list of class functions [f1,..,f_m] on H
# We output the mxn matrix with entries m_{ij} for which f_i = \sum_j m_{ij}ψj 
function rep_ring_decomposition_matrix(R, target_irr, class_functions)
    tbl = GAP.Globals.UnderlyingCharacterTable(target_irr[1])
    mults = GAP.Globals.MatScalarProducts(tbl, GapObj(target_irr), GapObj(class_functions))
    m = zero_matrix(R, length(class_functions), length(target_irr))
    for i in eachindex(class_functions), j in eachindex(target_irr)
        m[i, j] = R(mults[i][j])
    end
    return m
end

# Since the value of R_C at H is the free module on the irreducible characters of H, this will let us get res/tr as module maps


# res^K_H : R(K) -> R(H)
function rep_ring_restriction(R, M1, M2, irr1, irr2, H1)
    restricted = [GAP.Globals.RestrictedClassFunction(chi, H1) for chi in irr2]
    ModuleHomomorphism(M2, M1, rep_ring_decomposition_matrix(R, irr1, restricted))
end

# ind_H^K : R(H) -> R(K)
function rep_ring_transfer(R, M1, M2, irr1, irr2, H2)
    induced = [GAP.Globals.InducedClassFunction(chi, H2) for chi in irr1]
    ModuleHomomorphism(M1, M2, rep_ring_decomposition_matrix(R, irr2, induced))
end

# c_g : R(H) -> R(gHg^-1), sending chi to the character y ↦ chi(g^-1 y g)
function rep_ring_conjugation(R, mc, irreducibles, gi, Hi, values)
    g = mc.generators[gi]
    gHi = mc.generator_left_conjugation_matrix[gi, Hi]
    tbl = GAP.Globals.UnderlyingCharacterTable(irreducibles[gHi][1])
    # GAP writes y^g for g^-1*y*g
    reps = [GAP.Globals.Representative(C) for C in GAP.Globals.ConjugacyClasses(tbl)]
    conjugated = [
        GAP.Globals.ClassFunction(tbl, GapObj([(y^g)^chi for y in reps]))
        for chi in irreducibles[Hi]
    ]
    ModuleIsomorphism(values[Hi], values[gHi], rep_ring_decomposition_matrix(R, irreducibles[gHi], conjugated))
end

function _rep_ring_irreducibles(mc::MackeyContext)
    return [collect(GAP.Globals.Irr(GAP.Globals.CharacterTable(H))) for H in mc.subgroups]
end

"""
    representation_ring_mackey_functor(mc::MackeyContext, R::Ring = ZZ) -> MackeyFunctor

Return the complex representation ring Mackey functor for the group specified by `mc`, with coefficients in `R`.

The value at a subgroup ``H`` is ``R_{\\mathbb{C}}(H) \\otimes_{\\mathbb{Z}} R``, the free ``R``-module on the irreducible complex characters of ``H``. Restrictions are restrictions of characters, transfers are induced characters, and conjugation by ``g`` sends a character ``\\chi`` of ``H`` to the character ``y \\mapsto \\chi(g^{-1}yg)`` of ``gHg^{-1}``.
"""
function representation_ring_mackey_functor(mc::MackeyContext, R::Ring=ZZ)
    irreducibles = _rep_ring_irreducibles(mc)
    values = [free_module(R, length(irr)) for irr in irreducibles]

    cover_transfers = Generic.ModuleHomomorphism[
        rep_ring_transfer(R, values[i], values[j], irreducibles[i], irreducibles[j], mc.subgroups[j])
        for (i, j) in mc.covers
    ]
    cover_restrictions = Generic.ModuleHomomorphism[
        rep_ring_restriction(R, values[i], values[j], irreducibles[i], irreducibles[j], mc.subgroups[i])
        for (i, j) in mc.covers
    ]
    generator_conjugations = Generic.ModuleIsomorphism[
        rep_ring_conjugation(R, mc, irreducibles, gi, Hi, values)
        for gi in eachindex(mc.generators), Hi in eachindex(mc.subgroups)
    ]

    MackeyFunctor(mc, values, cover_restrictions, cover_transfers, generator_conjugations)
end

"""
    representation_ring_mackey_functor(G, R::Ring = ZZ) -> MackeyFunctor

Return the complex representation ring Mackey functor for the group `G`, with coefficients in `R`.
"""
representation_ring_mackey_functor(G::GapObj, R::Ring=ZZ) = representation_ring_mackey_functor(MackeyContext(G), R)
