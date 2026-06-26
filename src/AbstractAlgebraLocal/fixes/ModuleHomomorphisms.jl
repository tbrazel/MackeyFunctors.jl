# Checks if two module maps are identical (mathematically)
function map_eq(f::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism), g::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism))::Bool
    domain(f) === domain(g) && codomain(f) === codomain(g) && all(x -> f(x) == g(x), gens(domain(f)))
end

function _preserves_domain_relations(f::Generic.ModuleHomomorphism)::Bool
    # An FPModule is presented by generators modulo the rows returned by
    # `relations(domain(f))`.  A matrix gives a genuine map out of the quotient
    # only when each source relation is sent to zero in the target quotient.
    #
    # This check matters for isomorphism testing because AbstractAlgebra's
    # `ModuleIsomorphism` constructor solves for a matrix inverse in a larger
    # presentation matrix.  That inverse matrix can fail to respect the source
    # relations of the proposed inverse.  For example, the reduction map
    # Z/4 -> Z/2 admits the matrix [1] as a one-sided solution, but the candidate
    # inverse Z/2 -> Z/4 sends the relation 2 = 0 in Z/2 to 2 != 0 in Z/4.
    for relation in relations(domain(f))
        image_relation = relation * matrix(f)
        codomain(f)(image_relation) == zero(codomain(f)) || return false
    end
    return true
end

function direct_sum_homomorphism(
    source,
    target,
    p::Generic.ModuleIsomorphism,
    q::Generic.ModuleIsomorphism,
)
    return direct_sum_homomorphism(source, target, as_homomorphism(p), as_homomorphism(q))
end

"""
    direct_sum_homomorphism(source, target, p::Generic.ModuleHomomorphism, q::Generic.ModuleHomomorphism)

Return the direct sum of `p` and `q`. `source` and `target` are the domain and codomain of the direct sum morphism.
"""
function direct_sum_homomorphism(
    source,
    target,
    p::Generic.ModuleHomomorphism,
    q::Generic.ModuleHomomorphism,
)
    R = base_ring(source)
    all(
        mod -> base_ring(mod) == R,
        (target, domain(p), codomain(p), domain(q), codomain(q)),
    ) || throw(ArgumentError("All direct-sum homomorphism modules must have the same base ring."))

    source_generators = ngens(domain(p)) + ngens(domain(q))
    target_generators = ngens(codomain(p)) + ngens(codomain(q))
    ngens(source) == source_generators ||
        throw(ArgumentError("The source does not have the expected direct-sum presentation."))
    ngens(target) == target_generators ||
        throw(ArgumentError("The target does not have the expected direct-sum presentation."))

    block_matrix = zero_matrix(R, ngens(source), ngens(target))
    _copy_matrix_block!(block_matrix, matrix(p), 0, 0)
    _copy_matrix_block!(
        block_matrix,
        matrix(q),
        ngens(domain(p)),
        ngens(codomain(p)),
    )

    return ModuleHomomorphism(source, target, block_matrix)
end
