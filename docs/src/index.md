# MackeyFunctors.jl

`MackeyFunctors.jl` provides tools for constructing Mackey functors for arbitrary finite groups ``G``. Some supported functionality includes:
- a custom type for [Mackey functors](types/MackeyFunctors.md), which is sped up by precompiling some [information about the desired group](types/MackeyContexts.md)
- various [constructors](manual/constructors.md) including fixed point, free, Burnside, and shifts of arbitrary Mackey functors
- direct sums of Mackey functors

Planned future functionality includes:
- box products, internal and external homs
- free resolutions
- ext and tor computations
- support for cohomological Mackey functors

This project was first developed by a handful of people including a subset of the developers of the [`CpMackeyFunctors.m2`](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/CpMackeyFunctors/html/index.html) package for Macaulay2. Version 0.0 was developed at the [ICERM Machine Computation in Homotopy Theory](https://icerm.brown.edu/program/topical_workshop/tw-26-mch) conference in summer of 2026. We are grateful to ICERM and the organizers of this conference for the opportunity to work on this package.

## Benchmarks

```@eval
import Markdown
import BenchmarkTools
using GAP, AbstractAlgebra, MackeyFunctors

C2 = GAP.Globals.CyclicGroup(2)
S3 = GAP.Globals.SymmetricGroup(3)
S4 = GAP.Globals.SymmetricGroup(4)

C2_context = MackeyContext(C2)
S3_context = MackeyContext(S3)
S4_context = MackeyContext(S4)

benchmarks = [
    ("C2", "context", "`MackeyContext(C2)`",
        BenchmarkTools.@benchmarkable MackeyContext($C2) samples=5 seconds=0.5 evals=1),
    ("C2", "constant", "`constant_mackey_functor(C2_context, ZZ)`",
        BenchmarkTools.@benchmarkable constant_mackey_functor($C2_context, $ZZ) samples=5 seconds=0.5 evals=1),
    ("C2", "Burnside", "`burnside_mackey_functor(C2_context)`",
        BenchmarkTools.@benchmarkable burnside_mackey_functor($C2_context) samples=5 seconds=0.5 evals=1),
    ("S3", "context", "`MackeyContext(S3)`",
        BenchmarkTools.@benchmarkable MackeyContext($S3) samples=5 seconds=0.5 evals=1),
    ("S3", "constant", "`constant_mackey_functor(S3_context, ZZ)`",
        BenchmarkTools.@benchmarkable constant_mackey_functor($S3_context, $ZZ) samples=5 seconds=0.5 evals=1),
    ("S3", "Burnside", "`burnside_mackey_functor(S3_context)`",
        BenchmarkTools.@benchmarkable burnside_mackey_functor($S3_context) samples=5 seconds=0.5 evals=1),
    ("S4", "context", "`MackeyContext(S4)`",
        BenchmarkTools.@benchmarkable MackeyContext($S4) samples=5 seconds=0.5 evals=1),
    ("S4", "constant", "`constant_mackey_functor(S4_context, ZZ)`",
        BenchmarkTools.@benchmarkable constant_mackey_functor($S4_context, $ZZ) samples=5 seconds=0.5 evals=1),
    ("S4", "Burnside", "`burnside_mackey_functor(S4_context)`",
        BenchmarkTools.@benchmarkable burnside_mackey_functor($S4_context) samples=5 seconds=0.5 evals=1),
]

rows = [
    "The Mackey functor constructor rows use a precomputed `MackeyContext`; context construction is benchmarked separately.",
    "",
    "| Group | Operation | Method | Minimum time | Median time | Memory | Allocations |",
    "|:--|:--|:--|--:|--:|--:|--:|",
]

for (group, operation, label, b) in benchmarks
    trial = BenchmarkTools.run(b)
    min_trial = minimum(trial)
    med_trial = BenchmarkTools.median(trial)

    push!(rows,
        "| $group | $operation | $label | $(BenchmarkTools.prettytime(min_trial.time)) | " *
        "$(BenchmarkTools.prettytime(med_trial.time)) | " *
        "$(BenchmarkTools.prettymemory(min_trial.memory)) | $(min_trial.allocs) |"
    )
end

Markdown.parse(join(rows, "\n"))
```
