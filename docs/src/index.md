# MackeyFunctors.jl

`MackeyFunctors.jl` provides tools for constructing and validating
lattice-backed Mackey functor data for finite GAP groups.

The package represents values as finitely generated abelian groups using
invariant factors. A value vector `[0, 2, 4]` represents `Z x Z/2 x Z/4`.
Restriction, transfer, and conjugation maps are represented by integer
matrices between these presentations.

## Installation

From the package checkout:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Documentation Build

To build these docs locally:

```julia
using Pkg
Pkg.activate("docs")
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()
```

Then run:

```sh
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`.
