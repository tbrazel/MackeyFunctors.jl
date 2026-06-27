function _induced_submodule_homomorphism(
    f::Generic.ModuleHomomorphism,
    source_submodule,
    target_submodule,
)::Generic.ModuleHomomorphism
    source, = source_submodule
    target, = target_submodule
    if ngens(source) == 0
        return zero_homomorphism(source, target)
    end

    return ModuleHomomorphism(
        source,
        target,
        submodules_matrix(f, source_submodule, target_submodule),
    )
end

function _induced_submodule_isomorphism(
    f::Generic.ModuleHomomorphism,
    source_submodule,
    target_submodule,
)::Generic.ModuleIsomorphism
    homomorphism = _induced_submodule_homomorphism(
        f,
        source_submodule,
        target_submodule,
    )

    return ModuleIsomorphism(
        domain(homomorphism),
        codomain(homomorphism),
        matrix(homomorphism),
    )
end

"""
    kernel(f::MackeyFunctorHomomorphism) -> (MackeyFunctor, MackeyFunctorHomomorphism)

Return the kernel Mackey functor of `f`, together with its inclusion into
`f.domain`.

The construction is componentwise: the value at `H` is the module kernel of
`f.components[H]`, and the Mackey structure maps are induced from the domain
Mackey functor.
"""
function kernel(f::MackeyFunctorHomomorphism)
    context = f.context
    kernel_data = [kernel(component) for component in f.components]
    values = AbstractAlgebra.FPModule[first(data) for data in kernel_data]

    cover_restrictions = Generic.ModuleHomomorphism[
        _induced_submodule_homomorphism(
            f.domain.cover_restrictions[cover_index],
            kernel_data[j],
            kernel_data[i],
        )
        for (cover_index, (i, j)) in enumerate(context.covers)
    ]
    cover_transfers = Generic.ModuleHomomorphism[
        _induced_submodule_homomorphism(
            f.domain.cover_transfers[cover_index],
            kernel_data[i],
            kernel_data[j],
        )
        for (cover_index, (i, j)) in enumerate(context.covers)
    ]
    generator_conjugations = Generic.ModuleIsomorphism[
        _induced_submodule_isomorphism(
            as_homomorphism(f.domain.generator_conjugations[generator_index, H_index]),
            kernel_data[H_index],
            kernel_data[context.generator_left_conjugation_matrix[generator_index, H_index]],
        )
        for generator_index in eachindex(context.generators), H_index in eachindex(context.subgroups)
    ]

    kernel_mackey_functor = MackeyFunctor(
        context,
        values,
        cover_restrictions,
        cover_transfers,
        generator_conjugations,
    )
    inclusion = MackeyFunctorHomomorphism(
        kernel_mackey_functor,
        f.domain,
        Generic.ModuleHomomorphism[last(data) for data in kernel_data],
    )

    return kernel_mackey_functor, inclusion
end
