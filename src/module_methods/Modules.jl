# Checks if two module maps are identical (mathematically)
function map_eq(f::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism), g::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism))::Bool
    domain(f) === domain(g) && codomain(f) === codomain(g) && all(x -> f(x) == g(x), gens(domain(f)))
end

function as_homomorphism(f::Generic.ModuleIsomorphism)::Generic.ModuleHomomorphism
    return ModuleHomomorphism(domain(f), codomain(f), matrix(f))
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
    ModuleHomomorphism(
        M, N, zero_matrix(
            base_ring(M), ngens(M), ngens(N)))
end

function is_zero_module_homomorphism(phi::Generic.ModuleHomomorphism)::Bool
    return all(x -> phi(x) == zero(codomain(phi)), gens(domain(phi)))
end

function is_identity_module_homomorphism(phi::AbstractAlgebra.Map(AbstractAlgebra.FPModuleHomomorphism))::Bool
    return domain(phi) === codomain(phi) && all(phi(x) == x for x in gens(domain(phi)))
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
