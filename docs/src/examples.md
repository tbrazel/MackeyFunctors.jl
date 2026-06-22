# Examples

This example constructs the Burnside Mackey functor for the cyclic group `C2`.

```julia
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

M = MackeyFunctor(C2, values, restrictions, transfers, conjugations)
```

For interactive visualization data, use:

```julia
visualizer_data(M)
visualizer_json(M)
```
