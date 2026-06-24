using GAP
using AbstractAlgebra
using MackeyFunctors

C2 = GAP.Globals.CyclicGroup(2)



function prime_p(n::Int)
    n <= 1 && return false
    n == 2 && return true
    n % 2 == 0 && return false

    for i in 3:2:isqrt(n)
        if n % i == 0
            return false
        end
    end
    return true
end

struct makeUnderlyingFreeMackeyFunctor
    p       :: Int
    context :: MackeyContext
    modules :: Vector{AbstractAlgebra.FPModule}
    res     :: Generic.ModuleHomomorphism
    tr      :: Generic.ModuleHomomorphism
    conj    :: Generic.ModuleHomomorphism

    function makeUnderlyingFreeMackeyFunctor(p::Int)
        @assert prime_p(p) "p must be prime"

        Cp      = GAP.Globals.CyclicGroup(p)
        context = MackeyContext(Cp)

        Z  = free_module(ZZ, 1)
        Zp = free_module(ZZ, p)


        res = hom(Z, Zp, matrix(ZZ, 1, p, ones(Int, p)))   


        tr  = hom(Zp, Z,  matrix(ZZ, p, 1, ones(Int, p)))  

        P = zeros(Int, p, p)
        for i in 1:p
            P[i, mod(i, p) + 1] = 1
        end
        γ = hom(Zp, Zp, matrix(ZZ, p, p, vec(P)))

        new(p, context, [Z, Zp], res, tr, γ)
    end          
end              





