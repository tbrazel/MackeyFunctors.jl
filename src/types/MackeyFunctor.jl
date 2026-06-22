using MackeyContext

struct MackeyFunctor
    context::MackeyContext
    cover_restrictions::Vector{ModuleHomomorphism}
    cover_transfers::Vector{ModuleHomomorphism}
    generator_conjugations::Array{ModuleHomomorphism,2}
    generator_inverse_conjugations::Array{ModuleHomomorphism,2}
end