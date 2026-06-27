"""
    MackeyFunctorHomMackeyFunctor(M::MackeyFunctor, N::MackeyFunctor)
    InternalHom(M::MackeyFunctor, N::MackeyFunctor)

Represent the internal hom in the category of Mackey functors.

Mathematically, when `M` and `N` are `G`-Mackey functors, the internal hom is
the `G`-Mackey functor ``[M,N]`` whose values are

```math
[M,N](G/H) = \\operatorname{Hom}(M_{G/H}, N).
```

Here ``M_{G/H}`` represents the ``G/H``-shift of `M`. When `H` is a subgroup
of `K`, the map ``G/H \\to G/K`` induces maps ``M_{G/H} \\leftrightarrow
M_{G/K}``; these yield the restriction and transfer maps of ``[M,N]``.
Likewise for conjugation.
"""
function MackeyFunctorHomMackeyFunctor(
    domain_mf::MackeyFunctor,
    codomain_mf::MackeyFunctor;
    verify::Bool=true,
)::MackeyFunctor
    domain_mf.context == codomain_mf.context ||
        throw(ArgumentError("Mackey functors must have the same Mackey context."))
    coefficient_ring(domain_mf) == coefficient_ring(codomain_mf) ||
        throw(ArgumentError("Mackey functors must be defined over the same coefficient ring."))

    context = domain_mf.context

    # The value of the internal Hom at G/H is the external Hom module
    #
    #     [M,N](G/H) = Hom(M_H, N),
    #
    # where M_H is the G/H-shift of M.  We construct all shifted functors first
    # so the shift maps used below have domains and codomains identical to the
    # Hom modules' stored Mackey functor objects.
    shifted_domains = MackeyFunctor[
        shift(domain_mf, H_index; verify=verify)
        for H_index in eachindex(context.subgroups)
    ]
    hom_modules = MackeyFunctorHomModule[
        Hom(shifted_domains[H_index], codomain_mf)
        for H_index in eachindex(context.subgroups)
    ]
    values = AbstractAlgebra.FPModule[
        underlying_module(H)
        for H in hom_modules
    ]

    cover_restrictions = Generic.ModuleHomomorphism[]
    cover_transfers = Generic.ModuleHomomorphism[]
    for (H_index, K_index) in context.covers
        # A cover H <= K gives a projection of orbits G/H -> G/K.
        #
        # Covariantly, this gives shift_transfer(M,H,K): M_H -> M_K.  Since
        # Hom is contravariant in its first variable, precomposition with this
        # map is the internal-Hom restriction
        #
        #     Hom(M_K,N) -> Hom(M_H,N).
        push!(
            cover_restrictions,
            precomposition_map(
                hom_modules[K_index],
                hom_modules[H_index],
                shift_transfer(domain_mf, H_index, K_index; verify=verify),
            ),
        )

        # The opposite span gives shift_restriction(M,H,K): M_K -> M_H.
        # Precomposition with that map goes in the other direction and is the
        # internal-Hom transfer
        #
        #     Hom(M_H,N) -> Hom(M_K,N).
        push!(
            cover_transfers,
            precomposition_map(
                hom_modules[H_index],
                hom_modules[K_index],
                shift_restriction(domain_mf, H_index, K_index; verify=verify),
            ),
        )
    end

    generator_conjugations = Matrix{Generic.ModuleIsomorphism}(
        undef,
        length(context.generators),
        length(context.subgroups),
    )
    for generator_index in eachindex(context.generators), H_index in eachindex(context.subgroups)
        g = context.generators[generator_index]
        target_H_index =
            context.generator_left_conjugation_matrix[generator_index, H_index]

        # The generator conjugation of a Mackey functor goes
        #
        #     [M,N](H) -> [M,N](gHg^-1).
        #
        # An element on the left is a map M_H -> N.  To land in
        # Hom(M_{gHg^-1},N), precompose with the inverse shift-conjugation map
        #
        #     M_{gHg^-1} -> M_H,
        #
        # which is shift_conjugation(M, g^-1, gHg^-1).
        conjugation_map = precomposition_map(
            hom_modules[H_index],
            hom_modules[target_H_index],
            shift_conjugation(
                domain_mf,
                g^-1,
                target_H_index;
                verify=verify,
            ),
        )
        generator_conjugations[generator_index, H_index] = ModuleIsomorphism(
            values[H_index],
            values[target_H_index],
            matrix(conjugation_map),
        )
    end

    return MackeyFunctor(
        context,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations;
        verify=verify,
    )
end

function InternalHom(
    domain_mf::MackeyFunctor,
    codomain_mf::MackeyFunctor;
    verify::Bool=true,
)::MackeyFunctor
    return MackeyFunctorHomMackeyFunctor(
        domain_mf,
        codomain_mf;
        verify=verify,
    )
end
