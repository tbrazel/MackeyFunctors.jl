# work around a bug in AbstractAlgebra for non-free modules
for S in (Generic.ModuleHomomorphism, Generic.ModuleIsomorphism), T in (Generic.ModuleHomomorphism, Generic.ModuleIsomorphism)
    function Base.:(==)(f::S, g::T)
        (domain(f) === domain(g) && codomain(f) === codomain(g)) || return false
        if domain(f) isa Generic.FreeModule && codomain(f) isa Generic.FreeModule
            matrix(f) == matrix(g)
        else
            all(x -> f(x) == g(x), gens(domain(f)))
        end

    end
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
