as_homomorphism(f::Generic.ModuleHomomorphism) = f

function as_homomorphism(f::Generic.ModuleIsomorphism)::Generic.ModuleHomomorphism
    return ModuleHomomorphism(domain(f), codomain(f), matrix(f))
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

function _has_two_sided_inverse(
    f::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism),
    f_inv::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism),
)::Bool
    f = as_homomorphism(f)
    f_inv = as_homomorphism(f_inv)

    # First verify that both matrices define honest maps on the presented
    # quotient modules.  Once that is true, checking equality on the source
    # generators is enough for the two composite maps below.
    _preserves_domain_relations(f) || return false
    _preserves_domain_relations(f_inv) || return false

    return f * f_inv == identity_homomorphism(domain(f)) && f_inv * f == identity_homomorphism(codomain(f))
end

function is_isomorphism(f::Generic.ModuleIsomorphism)::Bool
    try
        return _has_two_sided_inverse(f, inv(f))
    catch
        return false
    end
end

function is_isomorphism(f::Generic.ModuleHomomorphism)
    try
        f_iso = ModuleIsomorphism(domain(f), codomain(f), matrix(f))
        return _has_two_sided_inverse(f, inv(f_iso))
    catch
        return false
    end
end

function is_isomorphism(
    f::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism),
)::Bool
    try
        return _has_two_sided_inverse(f, inv(f))
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

function is_zero_module_homomorphism(phi::Generic.ModuleHomomorphism)
    return all(iszero ∘ phi, gens(domain(phi)))
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
        target_matrix[row_offset + row, column_offset + column] =
            source_matrix[row, column]
    end

    return target_matrix
end
