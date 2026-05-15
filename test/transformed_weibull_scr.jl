# transformed_weibull_scr.jl
# This script is for investigation, not for automated testing.
# Plots the graph of the integrand φˢ(y) f(φ(y)) φ′(y) in
# the expression
#          ∞             1
#    μₛ =  ∫ xˢ f(x) dx = ∫ φˢ(y) f(φ(y)) φ′(y) dy,
#         0             -1
# where f denotes the weibull pdf with scale parameter 1
# and shape parameter k > 0 and s = 2n - 1.

# How to run this script: 
#     Open the folder for the local development version in VS Code:
#     File > Open Folder... > open the folder for the local 
#     development version of my Julia package CustomGaussQuadrature, 
#     which is in the folder:
#     \RESEARCH - NUMERICAL METHODS\QUADRATURE\Custom GAUSS\
#     CustomGaussQuadrature - Julia package\CustomGaussQuadrature\
#
#     Start the Julia REPL:
#     Ctrl+Shift+P  (to bring up Command Palette)>  Julia: Start REPL

using Pkg
# Activate the package root one level above this script.
# This makes `using CustomGaussQuadrature` load the local checkout,
# not any installed copy from the default environment.
Pkg.activate(joinpath(@__DIR__, ".."))
using CustomGaussQuadrature
pathof(CustomGaussQuadrature)
# Make the internal functions `ϕ_fn` and `ϕ′_fn` 
# available in this script.
using CustomGaussQuadrature: ϕ_fn, ϕ′_fn

using Plots

a = 0.0;
b = Inf;

"""
weibull_pdf = weibull_pdf_fn(T, x, k)

Returns the weibull pdf with scale parameter 1
and shape parameter k > 0, evaluated at x.
"""
function weibull_pdf_fn(T, x, k)
    Tk = parse(T, string(k))
    T0 = convert(T,0)
    T1 = convert(T, 1)
    @assert Tk > T0
    @assert x ≥ T0
    Tk * x^(Tk - T1) * exp(- x^Tk)
end

# Consider the test case of evaluation of
#   μₛ = ∫ xˢ f(x) dx 
# by numerical integration, for quite large s.
# We know that 
#   μₛ = Γ(1 + s/k).
# It is reasonable to consider s = 2n - 1,
# where n is the number of Gauss quadrature nodes.


"""
integrand_fn(T, s, y, a, b, k)

Returns the integrand φˢ(y) f(φ(y)) φ′(y) in
the expression
          ∞             1
    μₛ =  ∫ xˢ f(x) dx = ∫ φˢ(y) f(φ(y)) φ′(y) dy,
          0            -1
where f denotes the weibull pdf with scale parameter 1
and shape parameter k > 0.
"""
function integrand_fn(T, s, y, a, b, k)
    Tk = parse(T, string(k))
    T0 = convert(T,0)
    T1 = convert(T, 1)
    @assert Tk > T0
    @assert -T1 < y < T1
    x = ϕ_fn(T, y, a, b)
    (x^s) * weibull_pdf_fn(T, x, k) * ϕ′_fn(T, y, a, b)
end


"""
plot_integrand_fn(T, s, length, a, b, k)

Plots the graph of the integrand φˢ(y) f(φ(y)) φ′(y) in
the expression
          ∞             1
    μₛ =  ∫ xˢ f(x) dx = ∫ φˢ(y) f(φ(y)) φ′(y) dy,
         0             -1
where f denotes the weibull pdf with scale parameter 1
and shape parameter k > 0. 
"""
function plot_integrand_fn(T, s, length, a, b, k)
delta = 2.0 /length
y_grid = range(-1.0+delta, 1.0-delta, length=length);
integrand_grid = integrand_fn.(T, s, y_grid, a, b, k)
plot(y_grid, integrand_grid, 
title="f is weibull pdf with k=$k & s=$s",
titlefont=font(10), legend=false)
    xlims!(-1.0, 1.0)
    y_hi = maximum(integrand_grid)
    ylims!(0.0, y_hi)
    xlabel!("y")
    ylabel!("φˢ(y) f(φ(y)) φ′(y)")
end

T = BigFloat;
n = 10;
s = 2*n - 1
k = 2.1;
length = 10000;
plot_integrand_fn(T, s, length, a, b, k)






