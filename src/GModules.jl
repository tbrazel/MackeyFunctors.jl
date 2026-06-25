# functionality for G-modules

function map_extension(m, generator_maps, word)
    for (gi, k) in word
        gg = k >= 0 ? generator_maps[gi] : inv(generator_maps[gi])
        for _ in 1:abs(k)
            m = gg * m
        end
    end
    m
end

"""
    GModule(context::MackeyContext, M::AbstractAlgebra.FPModule, generator_actions::AbstractVector{<:Generic.ModuleIsomorphism}; verify::Bool = true)

Construct the equivariant module with underlying module `M`.
The action of the group `context.group` is given on the generators in `context.generators` by the isomorphisms in `generator_actions`.
"""
struct GModule
    context::MackeyContext
    M::AbstractAlgebra.FPModule
    generator_actions::Vector{<:Generic.ModuleIsomorphism}

    function GModule(context::MackeyContext, M::AbstractAlgebra.FPModule, generator_actions::AbstractVector{<:Generic.ModuleIsomorphism}; verify::Bool=true)
        if verify
            all(context.generator_relations) do word
                m = map_extension(identity_homomorphism(M), generator_actions, word)
                is_identity_module_homomorphism(m)
            end || throw(ArgumentError("The given matrices do not define a group action"))
        end
        new(context, M, generator_actions)
    end
end

function permutation_matrix(M, g)
    ginv = inv(g)
    m = zero_matrix(coefficient_ring(M), rank(M), rank(M))
    for i in 1:rank(M)
        m[i, i^ginv] = 1
    end
    ModuleIsomorphism(M, M, m)
end

"""
    permutation_module(context::MackeyContext, R::Ring = ZZ) -> GModule

Return the canonical `GModule` for the permutation group `context.group`.
"""
function permutation_module(context::MackeyContext, R::Ring=ZZ)
    GAP.Globals.IsPermGroup(context.group) || throw(ArgumentError("The given group is not a permutation group"))
    n = GAP.Globals.LargestMovedPoint(context.group)
    M = FreeModule(R, n)
    generator_actions = [permutation_matrix(M, g) for g in context.generators]
    GModule(context, M, generator_actions)
end
