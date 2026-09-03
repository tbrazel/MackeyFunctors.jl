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

Benchmarks are run locally (never in CI) via `julia --project=benchmark benchmark/benchmarks.jl`,
which regenerates the table below in place. See [`benchmark/benchmarks.jl`](https://github.com/tbrazel/MackeyFunctors.jl/blob/main/benchmark/benchmarks.jl).

```@raw html
<!-- BENCHMARK_TABLE_START -->
```
The Mackey functor constructor rows use a precomputed `MackeyContext`; context construction is benchmarked separately.

*Last updated 2026-09-03, Julia 1.12.6, arm64-apple-darwin24.0.0.*

| Group | Operation | Method | Minimum time | Median time | Memory | Allocations |
|:--|:--|:--|--:|--:|--:|--:|
| C2 | context | `MackeyContext(C2)` | 408.625 μs | 427.583 μs | 385.73 KiB | 4684 |
| C2 | constant | `constant_mackey_functor(C2_context, ZZ)` | 44.292 μs | 45.416 μs | 62.85 KiB | 1865 |
| C2 | Burnside | `burnside_mackey_functor(C2_context)` | 131.417 μs | 138.833 μs | 130.55 KiB | 3584 |
| S3 | context | `MackeyContext(S3)` | 3.557 ms | 3.862 ms | 2.93 MiB | 51277 |
| S3 | constant | `constant_mackey_functor(S3_context, ZZ)` | 500.125 μs | 518.792 μs | 541.38 KiB | 16741 |
| S3 | Burnside | `burnside_mackey_functor(S3_context)` | 4.117 ms | 4.304 ms | 3.72 MiB | 97159 |
| S4 | context | `MackeyContext(S4)` | 56.208 ms | 57.006 ms | 36.47 MiB | 738076 |
| S4 | constant | `constant_mackey_functor(S4_context, ZZ)` | 4.860 ms | 4.896 ms | 4.70 MiB | 148505 |
| S4 | Burnside | `burnside_mackey_functor(S4_context)` | 124.976 ms | 131.259 ms | 119.98 MiB | 3925338 |
```@raw html
<!-- BENCHMARK_TABLE_END -->
```
