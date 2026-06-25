# Documentation Build

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

The generated local site is written to `docs/build/`.
