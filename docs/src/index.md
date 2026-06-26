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