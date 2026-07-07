# Manual Test-Script Environment

This environment is for exploratory and comparison scripts under `test/` that
use packages which are not runtime dependencies of `CustomGaussQuadrature`.

From the package root, initialize the environment with:

```julia
using Pkg
Pkg.activate("test/manual")
Pkg.develop(path=".")
Pkg.instantiate()
```

Then run a manual script from the Julia REPL with, for example:

```julia
include("test/try_QuadGK_pkg_scr.jl")
include("test/try_PolyChaos_pkg_scr.jl")
```
