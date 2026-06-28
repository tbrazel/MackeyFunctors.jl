# AbstractAlgebra.jl has some random functionality, but is not very extensive.
# Most importantly, it does have a way to generate random elmenents of a finitely
# generated ZZ-module.

using RandomExtensions

function RandomExtensions.make(M::AbstractAlgebra.FPModule, coeff_range::AbstractUnitRange{Int})
    RandomExtensions.Make(M, coeff_range)
end

function rand(rng::AbstractRNG, sp::RandomExtensions.SamplerTrivial{<:RandomExtensions.Make2{<:AbstractAlgebra.FPModule,<:AbstractUnitRange{Int}}})
    M, coeff_range = sp[][1:end]
    xs = gens(M)
    R = base_ring(M)
    R_make = make(R, coeff_range)
    sum(rand(rng, R_make)*x for x in xs)
end

rand(rng::AbstractRNG, M::AbstractAlgebra.FPModule, coeff_range::AbstractUnitRange{Int}) =
    rand(rng, make(M, coeff_range))

rand(M::AbstractAlgebra.FPModule, coeff_range) = rand(Random.default_rng(), M, coeff_range)