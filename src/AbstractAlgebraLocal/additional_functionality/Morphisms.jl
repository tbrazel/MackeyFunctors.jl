as_homomorphism(f::Generic.ModuleHomomorphism) = f

function as_homomorphism(f::Generic.ModuleIsomorphism)::Generic.ModuleHomomorphism
    return ModuleHomomorphism(domain(f), codomain(f), matrix(f))
end

is_invertible(f) = AbstractAlgebra.is_invertible(f)

is_invertible(f::Generic.ModuleIsomorphism) = true

function is_invertible(f::Generic.ModuleHomomorphism)
    try
        ModuleIsomorphism(domain(f), codomain(f), matrix(f))
        return true
    catch
        return false
    end
end

# Given a module M, returns its identity as a type Generic.ModuleIsomorphism
function identity_isomorphism(M::AbstractAlgebra.FPModule)::Generic.ModuleIsomorphism
    ModuleIsomorphism(M, M, identity_matrix(base_ring(M), ngens(M)))
end

# Given a module M, returns its identity as a type Generic.ModuleHomomorphism
function identity_homomorphism(M::AbstractAlgebra.FPModule)::Generic.ModuleHomomorphism
    as_homomorphism(identity_isomorphism(M))
end

# Given modules M and N, returns the zero map from M to N
function zero_homomorphism(M::AbstractAlgebra.FPModule, N::AbstractAlgebra.FPModule)::Generic.ModuleHomomorphism
    ModuleHomomorphism(M, N, zero_matrix(base_ring(M), ngens(M), ngens(N)))
end

function Base.iszero(f::Generic.ModuleHomomorphism)
    if codomain(f) isa Generic.FreeModule
        iszero(matrix(f))
    else
        all(iszero ∘ f, gens(domain(f)))
    end
end

function is_identity_module_homomorphism(phi::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism))
    return domain(phi) === codomain(phi) && all(x -> phi(x) == x, gens(domain(phi)))
end

"""
    $(@__MODULE__).submodules_matrix([f], (M1, f1), (M2, f2))

Let `f1` and `f2` be morphisms with domains `M1` and `M2`, where `f2` is assumed to be injective.
Also, let `f` be a morphism from `codomain(M1)` to `codomain(M2)`. This function returns a morphism
`g` from `M1` to `M2` such that `f2(g(m))` equals `f(f1(m))` for all `m` in `M1`.

If `f` is not given, it is assumed to be the identity map.
"""
submodules_matrix

submodules_matrix((M1, f1), (M2, f2)) = submodules_matrix(identity_map(codomain(f1)), (M1, f1), (M2, f2))

function submodules_matrix(f, (M1, f1), (M2, f2))
    @assert codomain(f1) === domain(f) && codomain(f) === codomain(f2)
    imgs = [preimage(f2, f(f1(v))) for v in generators(M1)]
    matrix([imgs[i][j] for i in eachindex(generators(M1)), j in eachindex(generators(M2))])
end

function _copy_matrix_block!(target_matrix, source_matrix, row_offset::Int, column_offset::Int)
    for row in 1:nrows(source_matrix), column in 1:ncols(source_matrix)
        target_matrix[row_offset+row, column_offset+column] =
            source_matrix[row, column]
    end

    return target_matrix
end

function _check_module_homomorphism_matrix(
    maps::AbstractMatrix{<:Generic.ModuleHomomorphism},
)
    isempty(maps) &&
        throw(ArgumentError("Cannot build a homomorphism from an empty matrix."))

    row_indices = collect(axes(maps, 1))
    column_indices = collect(axes(maps, 2))

    for row in row_indices, column in column_indices
        if domain(maps[row, column]) !== domain(maps[row, first(column_indices)])
            throw(ArgumentError("Entries in each matrix row must have the same domain."))
        end
        if codomain(maps[row, column]) !== codomain(maps[first(row_indices), column])
            throw(ArgumentError("Entries in each matrix column must have the same codomain."))
        end
    end

    return nothing
end

"""
    block_homomorphism(maps::AbstractMatrix{<:Generic.ModuleHomomorphism})

Build the homomorphism represented by a matrix of module homomorphisms. This
uses AbstractAlgebra's matrix convention: rows index domain summands and
columns index codomain summands. Thus `maps[i, j]` must be a map from the
`i`th domain summand to the `j`th codomain summand, and the result has type
`direct_sum(row domains...) -> direct_sum(column codomains...)`.
"""
function block_homomorphism(
    maps::AbstractMatrix{<:Generic.ModuleHomomorphism},
)::Generic.ModuleHomomorphism
    _check_module_homomorphism_matrix(maps)

    row_indices = collect(axes(maps, 1))
    column_indices = collect(axes(maps, 2))
    domain_summands = AbstractAlgebra.FPModule[
        domain(maps[row, first(column_indices)])
        for row in row_indices
    ]
    codomain_summands = AbstractAlgebra.FPModule[
        codomain(maps[first(row_indices), column])
        for column in column_indices
    ]

    domain_sum, _, domain_projections = direct_sum(domain_summands)
    codomain_sum, codomain_injections, = direct_sum(codomain_summands)

    return block_homomorphism(
        domain_sum,
        codomain_sum,
        domain_projections,
        codomain_injections,
        maps,
    )
end

function block_homomorphism(
    maps::AbstractVector{<:Generic.ModuleHomomorphism};
    orientation::Symbol=:column,
)::Generic.ModuleHomomorphism
    if orientation === :column
        return block_homomorphism(reshape(collect(maps), :, 1))
    elseif orientation === :row
        return block_homomorphism(reshape(collect(maps), 1, :))
    else
        throw(ArgumentError("orientation must be either :row or :column."))
    end
end

"""
    block_homomorphism(source, target, source_projections, target_injections, maps)

Build a block homomorphism using an existing direct-sum presentation. This is
the same row-domain/column-codomain convention as `block_homomorphism(maps)`,
but the source and target direct sums, along with their projections and
injections, are supplied by the caller.
"""
function block_homomorphism(
    source::AbstractAlgebra.FPModule,
    target::AbstractAlgebra.FPModule,
    source_projections::AbstractVector{<:Generic.ModuleHomomorphism},
    target_injections::AbstractVector{<:Generic.ModuleHomomorphism},
    maps::AbstractMatrix{<:Generic.ModuleHomomorphism},
)::Generic.ModuleHomomorphism
    _check_module_homomorphism_matrix(maps)

    row_indices = collect(axes(maps, 1))
    column_indices = collect(axes(maps, 2))
    length(source_projections) == length(row_indices) ||
        throw(ArgumentError("There must be one source projection for each matrix row."))
    length(target_injections) == length(column_indices) ||
        throw(ArgumentError("There must be one target injection for each matrix column."))

    result = zero_homomorphism(source, target)
    for (row_number, row) in enumerate(row_indices)
        source_projection = source_projections[row_number]
        domain(source_projection) === source ||
            throw(ArgumentError("A source projection has the wrong domain."))

        for (column_number, column) in enumerate(column_indices)
            target_injection = target_injections[column_number]
            codomain(target_injection) === target ||
                throw(ArgumentError("A target injection has the wrong codomain."))

            domain(maps[row, column]) === codomain(source_projection) ||
                throw(ArgumentError("A matrix entry has the wrong domain for its row."))
            codomain(maps[row, column]) === domain(target_injection) ||
                throw(ArgumentError("A matrix entry has the wrong codomain for its column."))

            result += source_projection * maps[row, column] * target_injection
        end
    end

    return result
end
