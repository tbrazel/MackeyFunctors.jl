using GAP
using AbstractAlgebra
using MackeyFunctors

C2 = GAP.Globals.CyclicGroup(2)
subs = Vector{GapObj}(GAP.Globals.AllSubgroups(C2))
e = subs[1]
g = GAP.Globals.GeneratorsOfGroup(C2)[1]

values = IdDict{GapObj, Vector{Int64}}(
    e => [0],
    C2 => [0, 0],
)

restrictions = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
    (e, C2) => [1 2],
)

transfers = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
    (e, C2) => [0; 1;;],
)

conjugations = IdDict{Tuple{GapObj, GapObj}, Matrix{Int64}}(
    (e, g) => [1;;],
)

M_Z = MackeyFunctor(C2, values, restrictions, transfers, conjugations)

F5 = GF(5)
M_F5 = MackeyFunctor(
    C2,
    values,
    restrictions,
    transfers,
    conjugations;
    coefficient_ring=F5,
)

Zx, x = polynomial_ring(ZZ, :x)
M_Zx = MackeyFunctor(Zx, C2, values, restrictions, transfers, conjugations)

println("Default coefficient ring: ", M_Z.coefficient_ring)
println("Finite-field coefficient ring: ", M_F5.coefficient_ring)
println("Polynomial coefficient ring: ", M_Zx.coefficient_ring)
println("Number of lattice subgroups: ", length(M_Z.subgroups))
println("Number of cover relations: ", length(M_Z.covers))
println("Number of stored restrictions: ", length(M_Z.restrictions))
println("Number of stored transfers: ", length(M_Z.transfers))
println("Number of stored conjugations: ", length(M_Z.conjugations))
