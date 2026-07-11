import Random
using RandomExtensions

"""
    rand([rng], make(context, coefficient_range, max_generators, max_relations))
    rand([rng], context, coefficient_range, max_generators, max_relations)

Construct a random finitely presented Mackey functor over `context` with
integer coefficients. The generators form a direct sum of between one and
`max_generators` free Mackey functors at randomly chosen subgroup levels.
Between zero and `max_relations` random elements of that sum are chosen, with
coordinates in `coefficient_range`; the result is the cokernel of the
universal map hitting those elements.

Pass an explicit random-number generator to `rand` for reproducible output.
"""
RandomExtensions.maketype(
    ::MackeyContext,
    ::AbstractUnitRange{Int},
    ::Int,
    ::Int,
) = MackeyFunctor

function _check_random_mackey_functor_parameters(
    coefficient_range::AbstractUnitRange{Int},
    max_generators::Int,
    max_relations::Int,
)
    isempty(coefficient_range) &&
        throw(ArgumentError("The coefficient range must be nonempty."))
    max_generators >= 1 ||
        throw(ArgumentError("max_generators must be positive."))
    max_relations >= 0 ||
        throw(ArgumentError("max_relations must be nonnegative."))
    return nothing
end

function Random.rand(
    rng::Random.AbstractRNG,
    context::MackeyContext,
    coefficient_range::AbstractUnitRange{Int},
    max_generators::Int,
    max_relations::Int,
)
    return rand(
        rng,
        make(
            context,
            coefficient_range,
            max_generators,
            max_relations,
        ),
    )
end

function Random.rand(
    context::MackeyContext,
    coefficient_range::AbstractUnitRange{Int},
    max_generators::Int,
    max_relations::Int,
)
    return rand(
        Random.default_rng(),
        context,
        coefficient_range,
        max_generators,
        max_relations,
    )
end

function _random_free_mackey_functor(
    rng::Random.AbstractRNG,
    context::MackeyContext,
    number_of_summands::Int,
)
    levels = rand(rng, eachindex(context.subgroups), number_of_summands)
    summands = MackeyFunctor[
        free_mackey_functor(context, level, ZZ)
        for level in levels
    ]
    return first(direct_sum(summands))
end

function _random_element(
    rng::Random.AbstractRNG,
    module_value::AbstractAlgebra.FPModule,
    coefficient_range::AbstractUnitRange{Int},
)
    generators = gens(module_value)
    isempty(generators) && return zero(module_value)

    return sum(
        (base_ring(module_value)(rand(rng, coefficient_range)) * generator
         for generator in generators);
        init=zero(module_value),
    )
end

function _random_presentation_map(
    rng::Random.AbstractRNG,
    generators::MackeyFunctor,
    coefficient_range::AbstractUnitRange{Int},
    number_of_relations::Int,
)
    if number_of_relations == 0
        relations = zero_mackey_functor(generators.context, coefficient_ring(generators))
        return _zero_mackey_functor_homomorphism(relations, generators)
    end

    relation_maps = MackeyFunctorHomomorphism[]
    for _ in 1:number_of_relations
        level = rand(rng, eachindex(generators.context.subgroups))
        element = _random_element(
            rng,
            value(generators, level),
            coefficient_range,
        )
        push!(relation_maps, universal_map(generators, level, element))
    end

    return block_homomorphism(relation_maps)
end

function Random.rand(
    rng::Random.AbstractRNG,
    sampler::Random.SamplerTrivial{
        <:RandomExtensions.Make4{
            MackeyFunctor,
            <:MackeyContext,
            <:AbstractUnitRange{Int},
            Int,
            Int,
        },
    },
)
    context, coefficient_range, max_generators, max_relations = sampler[][1:end]
    _check_random_mackey_functor_parameters(
        coefficient_range,
        max_generators,
        max_relations,
    )

    generators = _random_free_mackey_functor(
        rng,
        context,
        rand(rng, 1:max_generators),
    )
    presentation = _random_presentation_map(
        rng,
        generators,
        coefficient_range,
        rand(rng, 0:max_relations),
    )
    return first(cokernel(presentation))
end
