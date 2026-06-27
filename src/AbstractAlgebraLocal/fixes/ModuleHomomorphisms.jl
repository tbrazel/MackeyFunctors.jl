# work around a bug in AbstractAlgebra for non-free modules
for S in (Generic.ModuleHomomorphism, Generic.ModuleIsomorphism), T in (Generic.ModuleHomomorphism, Generic.ModuleIsomorphism)
    function Base.:(==)(f::S, g::T)
        (domain(f) === domain(g) && codomain(f) === codomain(g)) || return false
        if codomain(f) isa Generic.FreeModule
            matrix(f) == matrix(g)
        else
            all(x -> f(x) == g(x), gens(domain(f)))
        end
    end
end

# make sure that morphisms are well-defined
for mor in (:ModuleHomomorphism, :ModuleIsomorphism)
    @eval function $mor(M1::AbstractAlgebra.FPModule{T}, M2::AbstractAlgebra.FPModule{T}, m::MatElem{T}) where T<:RingElement
        all(rel -> iszero(M2(rel * m)), relations(M1)) ||
            throw(ArgumentError("The given assignments for generators do not preserve relations"))
        AbstractAlgebra.$mor(M1, M2, m)
    end

    @eval function $mor(M1::AbstractAlgebra.FPModule{T}, M2::AbstractAlgebra.FPModule{T}, v::Vector{S}) where {T<:RingElement,S<:AbstractAlgebra.FPModuleElem{T}}
        all(rel -> iszero(sum(splat(*), zip(rel, v))), relations(M1)) ||
            throw(ArgumentError("The given assignments for generators do not preserve relations"))
        AbstractAlgebra.$mor(M1, M2, v)
    end
end

"""
    direct_sum(fv::AbstractVector{Generic.ModuleIsomorphism{T}}; domain, codomain) where T <: RingElement -> Generic.ModuleIsomorphism{T}

Return the direct sum of the isomorphisms in `fv`.
The `domain` and `codomain` of the resulting isomorphism are constructed unless they are specified.
"""
function direct_sum(
    fv::AbstractVector{Generic.ModuleIsomorphism{T}};
    domain = first(direct_sum(map(domain, fv))),
    codomain = first(direct_sum(map(codomain, fv))),
) where T <: RingElement
    block_matrix = block_diagonal_matrix(map(matrix, fv))
    return ModuleIsomorphism(domain, codomain, block_matrix)  # we should be able to specify the inverse
end

"""
    direct_sum(fv::AbstractVector{Generic.ModuleHomomorphism{T}}; domain, codomain) where T <: RingElement -> Generic.ModuleHomomorphism{T}

Return the direct sum of the homomorphisms in `fv`.
The `domain` and `codomain` of the resulting homomorphism are constructed unless they are specified.
"""
function direct_sum(
    fv::AbstractVector{Generic.ModuleHomomorphism{T}};
    domain = first(direct_sum(map(domain, fv))),
    codomain = first(direct_sum(map(codomain, fv))),
) where T <: RingElement
    block_matrix = block_diagonal_matrix(map(matrix, fv))
    return ModuleHomomorphism(domain, codomain, block_matrix)
end

direct_sum(f::Union{Generic.ModuleHomomorphism{T}, Generic.ModuleIsomorphism{T}}...; kw...) where T <: RingElement = direct_sum(collect(f); kw...)
