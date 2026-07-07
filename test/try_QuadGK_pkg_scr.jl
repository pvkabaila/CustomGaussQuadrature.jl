using QuadGK
using SpecialFunctions

using CustomGaussQuadrature
using Printf

# Choose the arithmetic used by QuadGK.gauss. The weighted-rule constructor
# determines its working floating-point type from the finite endpoints passed
# to gauss. To use BigFloat arithmetic, comment out the Float64 line and
# uncomment the BigFloat line.
# T_gauss = Float64
T_gauss = BigFloat

# Used only when T_gauss is BigFloat.
BigFloat_precision_bits = 256

if T_gauss === BigFloat
    setprecision(BigFloat, BigFloat_precision_bits, base=2)
end

if T_gauss === BigFloat
    rtol_gauss_default = parse(BigFloat, "1.0e-30")
else
    rtol_gauss_default = 1.0e-10
end

# Parameters for the sine-power density comparison. Here m is a positive
# integer and n is the number of Gauss nodes.
m = 20
n = 8

# Both the inbuilt package Dates and the package Printf have a format command.
using Dates
date_time_now = now();
println("\n", Dates.Date(date_time_now),
    "   ", Dates.format(date_time_now, "HH:MM"))

include(joinpath(@__DIR__, "utilities_scr.jl"))

# pkg_dir is the package root directory, regardless of whether the
# package was loaded from a local path or from the registry.
pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
println("\n", "The package root is:")
println_wrap(pkg_dir)
print("\n")

# The aim is to compare like with like:
#
# 1. QuadGK.gauss constructs an n-point Gaussian quadrature rule in the
#    original x variable for the weighted integral
#
#        1
#        ∫ g(x) f(x) dx.
#        0
#
# 2. CustomGaussQuadrature using moment determinants constructs an n-point
#    Gaussian quadrature rule in the same x variable, using the raw moments
#    of the same weight function f.
#
# The density used for this comparison is the sine-power density
#
#        fₘ(x) = π sin(πx)ᵐ / B(1/2, (m + 1)/2),   0 ≤ x ≤ 1,
#
# and fₘ(x) = 0 otherwise. The parameter m is a positive integer. The
# moment formula used below is described in the LaTeX file
#
#   2026-07-07 Computing Moments for a Sine-Power Density - ChatGPT.tex
#
# in the "2026 - work on the next version" directory next to this package
# workspace. Since fₘ has finite support [0, 1], QuadGK.gauss can be called
# directly on that interval. No support transformation is used here.

function print_rule_comparison_fn(n, nodes_gauss, weights_gauss,
        nodes_other, weights_other, other_label)
    nodes_other_label = string("nodes_", other_label)
    weights_other_label = string("weights_", other_label)

    println("        nodes_gauss                    ", nodes_other_label)
    for i in 1:n
        @printf("%2d   %.16e   %.16e\n", i, nodes_gauss[i], nodes_other[i])
    end
    print("\n")

    println("       weights_gauss                  ", weights_other_label)
    for i in 1:n
        @printf("%2d   %.16e   %.16e\n", i, weights_gauss[i], weights_other[i])
    end
    print("\n")

    println("maximum(abs.(nodes_gauss .- ", nodes_other_label, ")) = ",
        maximum(abs.(nodes_gauss .- nodes_other)))

    println("maximum(abs.((weights_gauss .- ", weights_other_label, ") ./ ",
        weights_other_label, ")) = ",
        maximum(abs.((weights_gauss .- weights_other) ./ weights_other)), "\n")
end

function sine_power_pdf_fn(::Type{T}, x, m::Integer) where {T<:AbstractFloat}
    @assert m > 0
    T_x = convert(T, x)
    T_zero = zero(T)
    T_one = one(T)
    @assert T_zero <= T_x <= T_one

    T_m = convert(T, m)
    T_half = convert(T, 1) / convert(T, 2)
    normalizing_constant = convert(T, pi) / beta(T_half, (T_m + T_one) / convert(T, 2))
    return normalizing_constant * sin(convert(T, pi) * T_x)^m
end

function lnf_sine_power_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "sine-power density"
    m = materialize_integer_spec_fn(which_f[3])
    @assert m > 0

    T_x = convert(T, x)
    T_zero = zero(T)
    T_one = one(T)
    @assert T_zero <= T_x <= T_one

    if T_x == T_zero || T_x == T_one
        return -convert(T, Inf)
    end

    T_m = convert(T, m)
    T_half = T_one / convert(T, 2)
    T_pi = convert(T, pi)
    normalizing_constant = T_pi / beta(T_half, (T_m + T_one) / convert(T, 2))
    return log(normalizing_constant) + T_m * log(sin(T_pi * T_x))
end

# Moment computation follows the LaTeX note. Put Y = X - 1/2, so that
# sin(πX) = cos(πY) and -1/2 ≤ Y ≤ 1/2. The density of Y is symmetric,
# hence all odd centred moments vanish.
#
# Define the unnormalised centred integrals
#
#             1/2
#     D(n,m) = ∫ yⁿ cos(πy)ᵐ dy.
#            -1/2
#
# For m = 0:
#
#     D(2r,0) = 1 / ((2r + 1) 2²ʳ),    D(2r+1,0) = 0.
#
# For m = 1:
#
#     D(0,1) = 2/π,
#     D(2r,1) = 1/(π 2²ʳ⁻¹) - (2r)(2r-1)D(2r-2,1)/π².
#
# For m ≥ 2:
#
#     D(n,m) = ((m-1)/m)D(n,m-2) - n(n-1)D(n-2,m)/(m²π²).
#
function sine_power_centered_integral_D_fn(::Type{T}, n::Integer, m::Integer) where {T<:AbstractFloat}
    @assert n >= 0
    @assert m >= 0

    if isodd(n)
        return zero(T)
    end

    T_one = one(T)
    T_two = convert(T, 2)
    T_pi = convert(T, pi)

    if m == 0
        return T_one / (convert(T, n + 1) * T_two^n)
    elseif m == 1
        if n == 0
            return T_two / T_pi
        end
        term1 = T_one / (T_pi * T_two^(n - 1))
        term2 = convert(T, n * (n - 1)) * sine_power_centered_integral_D_fn(T, n - 2, 1) / T_pi^2
        return term1 - term2
    else
        T_m = convert(T, m)
        term1 = convert(T, m - 1) * sine_power_centered_integral_D_fn(T, n, m - 2) / T_m
        term2 = zero(T)
        if n >= 2
            term2 = convert(T, n * (n - 1)) *
                sine_power_centered_integral_D_fn(T, n - 2, m) / (T_m^2 * T_pi^2)
        end
        return term1 - term2
    end
end

# The centred moments are μ′ₙ = E[(X - 1/2)ⁿ] = D(n,m)/D(0,m).
function sine_power_centered_moment_fn(::Type{T}, n::Integer, m::Integer) where {T<:AbstractFloat}
    if isodd(n)
        return zero(T)
    end
    return sine_power_centered_integral_D_fn(T, n, m) /
        sine_power_centered_integral_D_fn(T, 0, m)
end

# Convert centred moments to raw moments using X = 1/2 + Y:
#
#     E(Xⁿ) = Σᵣ binomial(n, 2r) (1/2)ⁿ⁻²ʳ E(Y²ʳ),
#
# where the sum is over r = 0, …, floor(n/2). Odd centred moments are zero,
# so only even centred moments appear.
function sine_power_raw_moment_fn(::Type{T}, n::Integer, m::Integer) where {T<:AbstractFloat}
    @assert n >= 0

    T_one = one(T)
    T_two = convert(T, 2)
    moment = zero(T)
    for r in 0:div(n, 2)
        centered_order = 2 * r
        moment += convert(T, binomial(n, centered_order)) *
            (T_one / T_two)^(n - centered_order) *
            sine_power_centered_moment_fn(T, centered_order, m)
    end
    return moment
end

println("Sine-power density on [0, 1]")
println("f_m(x) = pi * sin(pi*x)^m / beta(1/2, (m+1)/2)")
println("Compare QuadGK.gauss with CustomGaussQuadrature using moment determinants.", "\n")

let
which_f = ["sine-power density", [0, 1], m]
a, b = which_f[2]

@assert m > 0
@assert n >= 1

rtol_gauss = rtol_gauss_default

moment_fn = function (::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "sine-power density"
    @assert s >= 0
    m = materialize_integer_spec_fn(which_f[3])
    return sine_power_raw_moment_fn(T, s, m)
end

weight_fn = x -> sine_power_pdf_fn(T_gauss, x, m)

T_a = materialize_scalar_spec_fn(T_gauss, a)
T_b = materialize_scalar_spec_fn(T_gauss, b)

@time nodes_gauss, weights_gauss =
    QuadGK.gauss(weight_fn, n, T_a, T_b; rtol=rtol_gauss)

@time nodes_momentdets_raw, weights_momentdets_raw =
    custom_gauss_quad_all_fn(moment_fn, which_f, n)
nodes_momentdets = T_gauss[convert(T_gauss, x) for x in nodes_momentdets_raw]
weights_momentdets = T_gauss[convert(T_gauss, w) for w in weights_momentdets_raw]

println("m = ", m)
println("n = ", n)
println("T_gauss = ", T_gauss)
println("precision(T_gauss) = ", precision(T_gauss), " bits")
println("rtol_gauss = ", rtol_gauss, "\n")

print_rule_comparison_fn(
    n,
    nodes_gauss,
    weights_gauss,
    nodes_momentdets,
    weights_momentdets,
    "momentdets",
)
end

println("Compare QuadGK.gauss with the CustomGaussQuadrature Stieltjes procedure.", "\n")

which_f = ["sine-power density", [0, 1], m]
a, b = which_f[2]

rtol_gauss = rtol_gauss_default

weight_fn = x -> sine_power_pdf_fn(T_gauss, x, m)

T_a = materialize_scalar_spec_fn(T_gauss, a)
T_b = materialize_scalar_spec_fn(T_gauss, b)

@time nodes_gauss, weights_gauss =
    QuadGK.gauss(weight_fn, n, T_a, T_b; rtol=rtol_gauss)

lnf_typed_fn = lnf_sine_power_pdf_fn
mu0 = 1

@time include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))
nodes_stieltjes = T_gauss[convert(T_gauss, x) for x in stieltjes_nodes]
weights_stieltjes = T_gauss[convert(T_gauss, w) for w in stieltjes_weights]

println("m = ", m)
println("n = ", n)
println("T_gauss = ", T_gauss)
println("precision(T_gauss) = ", precision(T_gauss), " bits")
println("rtol_gauss = ", rtol_gauss)
println("r = ", r)
println("stieltjes_nbits = ", stieltjes_nbits, "\n")

print_rule_comparison_fn(
    n,
    nodes_gauss,
    weights_gauss,
    nodes_stieltjes,
    weights_stieltjes,
    "stieltjes",
)
