# CustomGaussQuadrature

[![Build Status](https://github.com/pvkabaila/CustomGaussQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/pvkabaila/CustomGaussQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amaster)

A Julia package for computing Gauss quadrature nodes and weights for
non-classical weight functions.

Given a nonnegative weight function $f$, the $n$-point Gauss quadrature
approximation to $\int g(x)\,f(x)\,dx$ is

$$\sum_{i=1}^n \lambda_i \; g(\tau_i),$$

where $\tau_1, \dots, \tau_n$ are the nodes and $\lambda_1, \dots, \lambda_n$ are
the corresponding weights. This approximation is exact whenever $g$ is a polynomial
of degree $\le 2n - 1$. In statistical applications, $f$ is commonly a probability density function (pdf).

For weight functions $f$ that lead to nodes that are roots of classical orthogonal polynomials of a continuous variable, such as Legendre, Hermite and Laguerre polynomials, 
the nodes and weights are readily available. When this is not the case, 
the Gauss rule must be custom-made. This is done in two steps:

1. **Step 1** — Compute the relevant recursion coefficients.
2. **Step 2** — Use these coefficients to compute the nodes and weights. 
  

This package provides two methods for Step 1:
(a) **moment determinants** and (b) the **Stieltjes procedure**.
Step 2 is easily carried out using `GenericLinearAlgebra.jl`.

The precision used in Step 1 is chosen adaptively to maximize the number of accurate digits when the computed nodes and weights are ultimately converted to `Float64`, for use in further extensive `Float64` computations. 

For full mathematical background and detailed usage, see the
[User Manual](doc/User_Manual.md).

## Installation

```julia
using Pkg
Pkg.add("CustomGaussQuadrature")
```

Then load the package:

```julia
using CustomGaussQuadrature
```

## Identifying the weight function

The weight function $f$ is identified by the array `which_f` with components:\
(i) the name given to $f$ (a string),\
(ii) support interval of $f$ specified by a 2-vector of the endpoints of the support of the weight function and\
(iii) parameter vector (if any).

There are five **built-in** weight functions whose moment formulae and
log-weight functions are already implemented. These include the scaled chi pdf with positive integer parameter `m`, identified by <br />
`which_f = ["scaled chi pdf", [0, Inf], m]` 

For a **user-defined** weight function
(e.g. the Weibull pdf with scale parameter 1), the user must supply the required moment or log-weight
function, as shown in the examples below.

## Quick Start

### Method 1 — Moment determinants

When a closed-form formula for the $s$'th moment
$\mu_s = \int x^s f(x)\,dx$ is available, 
use `custom_gauss_quad_all_fn`. Its main inputs are:

- `moment_fn` — a function with inputs `T` (floating-point type), `which_f` and `s` (moment order), returning the $s$'th moment. 
- `which_f` — identifies the weight function (see above).
- `n` — the number of Gauss quadrature nodes.

**Built-in** weight function scaled chi pdf with parameter `m` = 160 & number of nodes `n` = 5:

```julia
m = 160
which_f = ["scaled chi pdf", [0, Inf], m]

n = 5
nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n)

nodes = convert(Vector{Float64}, nodes)
weights = convert(Vector{Float64}, weights)
```

**User-defined** weight function Weibull pdf with scale parameter 1 and shape parameter `γ` = 2 & number of nodes `n` = 9:

```julia
using SpecialFunctions

which_f = ["weibull pdf", [0, Inf], 2.0]

function moment_weibull_pdf_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    γ = which_f[3]
    @assert γ > 0
    T_γ = parse(T, string(γ))
    @assert s ≥ 0
    if s == 0
        return(convert(T, 1))
    end
    gamma(convert(T, 1) + convert(T, s) / T_γ)
end

n = 9
nodes, weights = custom_gauss_quad_all_fn(moment_weibull_pdf_fn, which_f, n)

nodes = convert(Vector{Float64}, nodes)
weights = convert(Vector{Float64}, weights)
```

### Method 2 — Stieltjes procedure

This is more widely applicable than the moment determinants method since it requires only the 
evaluation of `μ₀`, `μ₁` and $\log f(x)$.
For **built-in** weight functions, these are already
implemented in `stieltjes_lnf_stored_scr.jl`. For a **user-defined** weight function,
the user supplies `lnf_fn` (evaluating $\log f$) and a function for evaluating `μ₀` and `μ₁`. These scripts are invoked via `include()`.

**Built-in** weight function scaled chi pdf with parameter `m` = 160 & number of nodes `n` = 33:

```julia
m = 160
which_f = ["scaled chi pdf", [0, Inf], m]
n = 33

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
# Evaluation of stieltjes_nodes and stieltjes_weights:
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

nodes = convert(Vector{Float64}, stieltjes_nodes)
weights = convert(Vector{Float64}, stieltjes_weights)
```

**User-defined** weight function Weibull pdf with scale parameter 1 and shape parameter `γ` = 2 & number of nodes `n` = 10:


```julia
using SpecialFunctions

which_f = ["weibull pdf", [0, Inf], 2.0]
n = 10
T = BigFloat

function lnf_weibull_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    @assert x > convert(T, 0)
    γ = which_f[3]
    T_γ = parse(T, string(γ))
    @assert T_γ > convert(T, 0)
    log(T_γ) + (T_γ - convert(T, 1)) * log(x) - x^T_γ
end

function μ₀_μ₁_weibull_pdf_fn(::Type{T}, which_f) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    γ = which_f[3]
    T_γ = parse(T, string(γ))
    @assert T_γ > convert(T, 0)
    T_1 = convert(T, 1)
    [convert(T, 1) gamma(T_1 + T_1 / T_γ)]
end

lnf_fn = x -> lnf_weibull_pdf_fn(T, which_f, x)
μ₀, μ₁ = μ₀_μ₁_weibull_pdf_fn(T, which_f)

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
# Evaluation of stieltjes_nodes and stieltjes_weights:
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

nodes = convert(Vector{Float64}, stieltjes_nodes)
weights = convert(Vector{Float64}, stieltjes_weights)
```
