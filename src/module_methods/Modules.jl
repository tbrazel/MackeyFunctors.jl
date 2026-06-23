
# Given a module M, returns its identity as a type Generic.ModuleIsomorphism
function identity_isomorphism(M::AbstractAlgebra.FPModule)::Generic.ModuleIsomorphism
    ModuleIsomorphism(M, M, identity_matrix(base_ring(M), ngens(M)))
end

function zero_homomorphism(M::AbstractAlgebra.FPModule,N::AbstractAlgebra.FPModule)
    ModuleHomomorphism(M,N,zero_matrix(base_ring(M), ngens(N),ngens(M)))
end

function is_zero_module_homomorphism(phi::Generic.ModuleHomomorphism)
    return all(x -> phi(x) == zero(codomain(phi)), gens(domain(phi)))
end

function is_equal_module_homomorphism(phi::Generic.ModuleHomomorphism, psi::Generic.ModuleHomomorphism)
    return domain(phi) === domain(psi) && codomain(phi) === codomain(psi) && is_zero_module_homomorphism(phi - psi)
end

function is_identity_module_homomorphism(phi::Generic.ModuleIsomorphism)
    return domain(phi) === codomain(phi) && all(phi(x) == identity_isomorphism(domain(phi))(x) for x in gens(domain(phi)))
end

# Checks if two module maps are identical (mathematically)
function same_module_map(f, g)
    domain(f) === domain(g) || return false
    codomain(f) === codomain(g) || return false
    return all(x -> f(x) == g(x), gens(domain(f)))
end

function toHomomorphism(f::Generic.ModuleIsomorphism)
    return ModuleHomomorphism(domain(f),codomain(f),matrix(f))
end