# Workaround for AbstractAlgebra direct_sum relation rows on quotient modules.
#
# AbstractAlgebra 0.50.1 has correct element arithmetic for some direct sums of
# finitely presented modules, but the stored relation rows can be wrong when
# quotient summands are repeated.  Until that is fixed upstream, all package
# code that needs relation data from a direct sum should go through these
# helpers instead of AbstractAlgebra.direct_sum.

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

function _typed_direct_sum_summands(
    summands::AbstractVector{<:AbstractAlgebra.FPModule},
)
    isempty(summands) &&
        throw(ArgumentError("Cannot take the direct sum of no modules."))

    R = base_ring(first(summands))
    all(summand -> base_ring(summand) == R, summands) ||
        throw(ArgumentError("Direct-sum summands must have the same base ring."))

    T = elem_type(R)
    return AbstractAlgebra.FPModule{T}[summands...]
end

function _direct_sum_presentation(
    summands::Vector{<:AbstractAlgebra.FPModule{T}},
) where T <: RingElement
    isempty(summands) &&
        throw(ArgumentError("Cannot take the direct sum of no modules."))

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

    return direct_sum_module, projection_from_free, free_sum, offsets
end

function _direct_sum_module(
    summands::AbstractVector{<:AbstractAlgebra.FPModule},
)
    typed_summands = _typed_direct_sum_summands(summands)
    direct_sum_module, _, _, _ = _direct_sum_presentation(typed_summands)
    return direct_sum_module
end

function _direct_sum(summands::AbstractVector{<:AbstractAlgebra.FPModule})
    typed_summands = _typed_direct_sum_summands(summands)
    return _direct_sum(typed_summands)
end

function _direct_sum(
    summands::Vector{<:AbstractAlgebra.FPModule{T}},
) where T <: RingElement
    direct_sum_module, projection_from_free, free_sum, offsets =
        _direct_sum_presentation(summands)

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
