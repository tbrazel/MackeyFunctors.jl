# The component of the linearization map at H : A(H) -> R_C(H).
# The basis element of A(H) indexed by the H-conjugacy class of L <= H is the H-set H/L,
# which is sent to its permutation character C[H/L] = Ind_L^H(1).
function linearization_component(R, M1, M2, conj_classes, irr, H)
    permutation_characters = [
        GAP.Globals.PermutationCharacter(H, GAP.Globals.Representative(LL))
        for LL in conj_classes
    ]
    ModuleHomomorphism(M1, M2, rep_ring_decomposition_matrix(R, irr, permutation_characters))
end

"""
    linearization(mc::MackeyContext, R::Ring = ZZ) -> MackeyFunctorHomomorphism

Return the linearization homomorphism from the [`burnside_mackey_functor`](@ref) to the [`representation_ring_mackey_functor`](@ref) for the group specified by `mc`, with coefficients in `R`.

At level ``H``, this sends a finite ``H``-set ``X`` to its permutation representation ``\\mathbb{C}[X]``. On the basis of ``A(H)``, the transitive ``H``-set ``H/L`` is sent to the permutation character ``\\mathrm{Ind}_L^H(1)``.
"""
function linearization(mc::MackeyContext, R::Ring=ZZ)
    return linearization(burnside_mackey_functor(mc, R), representation_ring_mackey_functor(mc, R))
end

"""
    linearization(burnside::MackeyFunctor, rep_ring::MackeyFunctor) -> MackeyFunctorHomomorphism

Return the linearization homomorphism from `burnside` to `rep_ring`, which should be the outputs of [`burnside_mackey_functor`](@ref) and [`representation_ring_mackey_functor`](@ref) for the same Mackey context and coefficient ring.
"""
function linearization(burnside::MackeyFunctor, rep_ring::MackeyFunctor)
    mc = burnside.context
    mc == rep_ring.context || throw(ArgumentError("The Burnside and representation ring Mackey functors must have the same Mackey context."))
    R = coefficient_ring(burnside)
    R == coefficient_ring(rep_ring) || throw(ArgumentError("The Burnside and representation ring Mackey functors must have the same coefficient ring."))

    conj_classes = _burnside_conjugacy_classes(mc)
    irreducibles = _rep_ring_irreducibles(mc)
    for i in eachindex(mc.subgroups)
        rank(burnside.values[i]) == length(conj_classes[i]) ||
            throw(ArgumentError("The first argument must be a Burnside Mackey functor."))
        rank(rep_ring.values[i]) == length(irreducibles[i]) ||
            throw(ArgumentError("The second argument must be a representation ring Mackey functor."))
    end

    components = Generic.ModuleHomomorphism[
        linearization_component(R, burnside.values[i], rep_ring.values[i], conj_classes[i], irreducibles[i], H)
        for (i, H) in enumerate(mc.subgroups)
    ]

    MackeyFunctorHomomorphism(burnside, rep_ring, components)
end

"""
    linearization(G, R::Ring = ZZ) -> MackeyFunctorHomomorphism

Return the linearization homomorphism from the Burnside Mackey functor to the complex representation ring Mackey functor for the group `G`, with coefficients in `R`.
"""
linearization(G::GapObj, R::Ring=ZZ) = linearization(MackeyContext(G), R)
