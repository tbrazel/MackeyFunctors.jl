"""
    MackeyFunctorHomomorphism(ctx,D,C,comps)

todo
"""
struct MackeyFunctorHomomorphism
    context::MackeyContext
    domain::MackeyFunctor
    codomain::MackeyFunctor
    components::Vector{Generic.ModuleHomomorphism}

    function MackeyFunctorHomomorphism(
        domain_mf::MackeyFunctor,
        codomain_mf::MackeyFunctor,
        components::AbstractVector{<:Generic.ModuleHomomorphism},
    )
        domain_mf.context == codomain_mf.context || throw(ArgumentError("The domain and codomain Mackey functors must have the same Mackey context."))

        # Make sure the number of inputted component maps matches the number of subgroups
        ctx = codomain_mf.context
        length(components) == length(ctx.subgroups) || throw(ArgumentError("There must be one component map for each subgroup."))

        # Make sure the ith component map starts and ends at the value of our Mackey functors on the ith subgroup
        for i in eachindex(ctx.subgroups)
            component = components[i]
            domain(component) === domain_mf.values[i] || throw(ArgumentError("Component $i has the wrong domain."))
            codomain(component) === codomain_mf.values[i] || throw(ArgumentError("Component $i has the wrong codomain."))
        end

        # For each cover H<K, assert that the cover res/tr commute with the levels of our maps
        for (n, (i, j)) in enumerate(ctx.covers)
            #=
                         components[i]
                M_dom(H[i]) ----> M_codom(H[i])
                    |                   |
                tr  |                   |  tr
                    ↓                   ↓
                M_dom(H[j]) ----> M_codom(H[j])
                         components[j]
            =#
            components[i] * codomain_mf.cover_transfers[n] == domain_mf.cover_transfers[n] * components[j] ||
                throw(ArgumentError("Cover transfers don't commute with the component maps."))

            #=
                         components[i]
                M_dom(H[i]) ----> M_codom(H[i])
                    ↑                   ↑
                res |                   |  res
                    |                   |
                M_dom(H[j]) ----> M_codom(H[j])
                         components[j]
            =#
            domain_mf.cover_restrictions[n] * components[i] == components[j] * codomain_mf.cover_restrictions[n] ||
                throw(ArgumentError("Cover restrictions don't commute with the component maps."))
        end

        # Check generator conjugation commutes with component maps

        for generator_index in eachindex(ctx.generators), subgp_index in eachindex(ctx.subgroups)
            #=
                         component_H
                M_dom(H)      -->     M_codom(H)
                    |                   |
                c_g |                   |  c_g
                    ↓                   ↓
                M_dom(gHg^-1))--> M_codom(gHg^-1)
                         component_gHg^-1
            =#

            gHginvs_index = ctx.generator_left_conjugation_matrix[generator_index, subgp_index]

            (domain_mf.generator_conjugations[generator_index, subgp_index] * components[gHginvs_index] ==
             components[subgp_index] * codomain_mf.generator_conjugations[generator_index, subgp_index]) ||
                throw(ArgumentError("Generator conjugations don't commute with component maps"))
        end

        component_maps = Generic.ModuleHomomorphism[components...]

        return new(
            ctx,
            domain_mf,
            codomain_mf,
            component_maps
        )
    end
end

function id_homomorphism(mf::MackeyFunctor)::MackeyFunctorHomomorphism
    return MackeyFunctorHomomorphism(
        mf,
        mf,
        [identity_homomorphism(M) for M in mf.values]
    )
end

AbstractAlgebraLocal.is_invertible(f::MackeyFunctorHomomorphism) = all(is_invertible, f.components)

function AbstractAlgebra.compose(f::MackeyFunctorHomomorphism, g::MackeyFunctorHomomorphism)::MackeyFunctorHomomorphism
    f.codomain == g.domain || throw(ArgumentError("The codomain of the first homomorphism must equal the domain of the second homomorphism."))
    return MackeyFunctorHomomorphism(
        f.domain,
        g.codomain,
        [f.components[i] * g.components[i] for i in eachindex(f.components)]
    )
end