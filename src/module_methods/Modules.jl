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