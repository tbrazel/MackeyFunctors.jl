function _module_homomorphism_from_images(
    source::AbstractAlgebra.FPModule,
    target::AbstractAlgebra.FPModule,
    images::Vector,
)::Generic.ModuleHomomorphism
    if ngens(source) == 0
        return zero_homomorphism(source, target)
    end

    return ModuleHomomorphism(source, target, images)
end

function _induced_quotient_homomorphism(
    source_projection::Generic.ModuleHomomorphism,
    target_projection::Generic.ModuleHomomorphism,
    f::Generic.ModuleHomomorphism,
)::Generic.ModuleHomomorphism
    source = codomain(source_projection)
    target = codomain(target_projection)
    images = elem_type(target)[
        target_projection(f(preimage(source_projection, x)))
        for x in gens(source)
    ]

    return _module_homomorphism_from_images(source, target, images)
end

function _induced_quotient_isomorphism(
    source_projection::Generic.ModuleHomomorphism,
    target_projection::Generic.ModuleHomomorphism,
    f::Generic.ModuleHomomorphism,
)::Generic.ModuleIsomorphism
    homomorphism = _induced_quotient_homomorphism(
        source_projection,
        target_projection,
        f,
    )

    return ModuleIsomorphism(
        domain(homomorphism),
        codomain(homomorphism),
        matrix(homomorphism),
    )
end

"""
    cokernel(f::MackeyFunctorHomomorphism) -> (MackeyFunctor, MackeyFunctorHomomorphism)

Return the cokernel Mackey functor of `f`, together with the quotient map
from `f.codomain`.

The construction is componentwise: the value at `H` is the module quotient
`f.codomain(H) / image(f.components[H])`, and the Mackey structure maps are
induced from the codomain Mackey functor.
"""
function cokernel(f::MackeyFunctorHomomorphism)
    context = f.context
    image_data = [image(component) for component in f.components]
    quotient_data = [
        quo(f.codomain.values[i], first(image_data[i]))
        for i in eachindex(context.subgroups)
    ]
    values = AbstractAlgebra.FPModule[first(data) for data in quotient_data]
    quotient_maps = Generic.ModuleHomomorphism[last(data) for data in quotient_data]

    cover_restrictions = Generic.ModuleHomomorphism[
        _induced_quotient_homomorphism(
            quotient_maps[j],
            quotient_maps[i],
            f.codomain.cover_restrictions[cover_index],
        )
        for (cover_index, (i, j)) in enumerate(context.covers)
    ]
    cover_transfers = Generic.ModuleHomomorphism[
        _induced_quotient_homomorphism(
            quotient_maps[i],
            quotient_maps[j],
            f.codomain.cover_transfers[cover_index],
        )
        for (cover_index, (i, j)) in enumerate(context.covers)
    ]
    generator_conjugations = Generic.ModuleIsomorphism[
        _induced_quotient_isomorphism(
            quotient_maps[H_index],
            quotient_maps[context.generator_left_conjugation_matrix[generator_index, H_index]],
            as_homomorphism(f.codomain.generator_conjugations[generator_index, H_index]),
        )
        for generator_index in eachindex(context.generators), H_index in eachindex(context.subgroups)
    ]

    cokernel_mackey_functor = MackeyFunctor(
        context,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations,
    )
    projection = MackeyFunctorHomomorphism(
        f.codomain,
        cokernel_mackey_functor,
        quotient_maps,
    )

    return cokernel_mackey_functor, projection
end
