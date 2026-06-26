# Installation

When running the package for the first time, clone the repository and then navigate to the directory `MackeyFunctors.jl` locally on your machine. Run `julia` to open a Julia session, then run the following commands:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
Activate will activate the repository, then instantiate will download any and all dependencies you might need.

After this, any time you want to run the package, you should activate the package (this sets your context and gets everything ready to go), then load in the dependencies before you start coding. This looks like:
```julia
using Pkg
Pkg.activate(".")
using GAP, AbstractAlgebra, MackeyFunctors
# then run whatever you want here!
```
If you want to get started for the first time and see what this package can do, we recommend starting out with the [examples](examples.md).