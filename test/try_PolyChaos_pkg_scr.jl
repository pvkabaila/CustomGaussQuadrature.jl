using PolyChaos
using SpecialFunctions

using CustomGaussQuadrature
using Printf

# Both the inbuilt package Dates and the package Printf have a format command
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

# The aim is to compare Gauss rules, then the right PolyChaos object is not
# the raw output of quadgp itself.  quadgp gives a general-purpose quadrature
# rule.  To obtain PolyChaos's Gauss rule for the same weight, construct an
# OrthoPoly with addQuadrature=true; this uses quadgp internally for an
# arbitrary weight.

# According to
# Julia PACKAGES\PolyChaos\PolyChaos - SciML documentation\10. Functions
# the function 
# quadgp 
# computes the N-point quadrature rule for weight with support (lb, ub). 
# It does so by implementing the same transformation of the support interval
# to (-1,1) as that described on p.94 of 
# Gautschi (2004) Orthogonal Polynomials: Computation and Approximation
# The default N-point approximation to the integral is obtained by the 
# clenshaw_curtis rule. 

println("Example 1: weight is exp(-|x|ᵏ) for x ∈ (-∞, ∞), where the parameter k > 0")

# A let block gives this example its own local scope, so simple names
# such as which_f, n, k, nodes_polychaos, and weights_polychaos can be
# reused in later examples without accidental cross-example dependence.
let
    which_f = ["exp(-abs(x)^k)", [-Inf, Inf], "0.5"]
# All example-specific data for this weight is stored in which_f.  The
# Float64 value k below is derived from which_f only for the PolyChaos
# weight closure.
a, b = which_f[2]
k = materialize_scalar_spec_fn(Float64, which_f[3])
@assert k > 0

n = 5

println("Compare the PolyChaos Gauss rule with")
println("CustomGaussQuadrature using moment determinants.", "\n")

# For the weight w(x) = exp(-|x|^k), k > 0, the s-th moment is
#
#        ∞
#   μₛ = ∫ xˢ exp(-|x|ᵏ) dx.
#       -∞
#       
# Since w is even, the integrand is odd when s is odd, so μₛ = 0.
# For even s, symmetry gives
#
#         ∞  
#   μₛ = 2 ∫ xˢ exp(-xᵏ) dx.
#         0
#
# With u = xᵏ, x = u^(1/k), and
# dx = (1/k) u^(1/k - 1) du. Hence
#
#             ∞
#   μₛ = (2/k) ∫ u^((s + 1)/k - 1) exp(-u) du
#             0
#  
#       = (2/k) Γ((s + 1)/k).
#
moment_fn = function (::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "exp(-abs(x)^k)"
    @assert s ≥ 0
    T_k = materialize_scalar_spec_fn(T, which_f[3])
    @assert T_k > zero(T)
    if isodd(s)
        return zero(T)
    end
    T_s = convert(T, s)
    return convert(T, 2) * gamma((T_s + one(T)) / T_k) / T_k
end

# The four inputs to Measure are:
# (1) a name for the measure,
# (2) the weight function,
# (3) the support interval,
# (4) a symmetry indicator.
#
# The last input is true if the weight function is symmetric, i.e.
# w(-x) = w(x), and false otherwise.  Here the weight function is
# exp(-|x|^k), which is symmetric about 0, so the symmetry indicator is true.
#
# For the installed PolyChaos v1.2.0, the symmetry indicator is stored, then
# used mainly as a zero-detection shortcut for scalar products.
# In the source, Measure(..., symm) just saves the Boolean as m.symmetric;
# issymmetric(m) later returns that value. When OrthoPoly builds recurrence
# coefficients, however, the flag is not used: rm_compute(measure, ...) passes
# only m.w and m.dom into the quadrature/Stieltjes machinery. So it does not
# enforce symmetry, fold the interval, or explicitly set recurrence
# coefficients such as α_n = 0.
measure = Measure("ex1_meas", x -> exp(-(abs(x)^k)), (a, b), true)


# PolyChaos uses an all-purpose rule on the transformed interval [-1,1]
# to approximate the integrals needed by its recurrence-coefficient
# computation. By default that all-purpose rule is clenshaw_curtis.
# If we change it to fejer, we are changing the numerical machinery used
# to estimate the recurrence coefficients, not asking PolyChaos to return
# a Fejer rule at the end.
#
# The overall pipeline is:
# 1. Start from the measure, here w(x) = exp(-|x|^k) on (-Inf, Inf).
# 2. Use quadrature = clenshaw_curtis or quadrature = fejer inside
#    PolyChaos's general-purpose integration step to approximate the
#    integrals needed for recurrence-coefficient construction.
# 3. Obtain approximate recurrence coefficients (α, β).
# 4. Build the final Gauss quadrature rule from those recurrence
#    coefficients.
# 5. Store that final Gauss rule inside the OrthoPoly object.
# 6. Extract it with nw.
#
# Therefore, if we wrote 
#   opq = OrthoPoly(
#       "ex1_op",
#       n,
#       measure;
#       Nquad = Nquad,
#       quadrature = fejer,
#       addQuadrature = true,
#   )
# the final object would still contain a Gauss rule since 
# the recurrence coefficients used to build that Gauss rule had been
# computed using Fejer as the internal auxiliary quadrature.
#
# The keyword Nquad controls how many nodes are used by this internal
# all-purpose rule. Increasing Nquad may improve the approximation to the
# recurrence coefficients, and hence may improve the final Gauss rule.
# The keyword addQuadrature controls whether PolyChaos stores the final
# Gauss rule in the OrthoPoly object so that it can later be extracted by nw.
#
# The relevant constructor signature in PolyChaos is:
# OrthoPoly(
#    name::String, deg::Int, measure::AbstractMeasure;
#    Nrec = deg + 1,
#    Nquad = 10 * Nrec,
#    quadrature::Function = clenshaw_curtis,
#    discretization::Function = stieltjes,
#    addQuadrature::Bool = true
# )
# Here Nrec is not used only indirectly through the default value of Nquad.
# It is passed directly into the OrthoPoly machinery and tells PolyChaos how
# many recurrence coefficients (α, β) to compute and store. The default
# choice Nrec = deg + 1 is the smallest natural choice consistent with a
# polynomial family represented up to degree deg.
#
# Nquad is different: it controls the size of the internal all-purpose
# quadrature used while PolyChaos numerically computes those recurrence
# coefficients. The default Nquad = 10 * Nrec ties the auxiliary quadrature
# size to the requested amount of recurrence information, but the two
# keywords can be set separately.
#
# A practical limitation to keep in mind is that PolyChaos is mixed-precision
# rather than uniformly generic. Some downstream routines can accept general
# real-valued recurrence coefficients, but the standard arbitrary-weight path
# used by OrthoPoly(..., measure; ...) builds those recurrence coefficients
# through Float64-centered machinery such as quadgp/rm_compute. This is
# therefore not a high-precision BigFloat analogue of the corresponding
# CustomGaussQuadrature workflow.

# The function OrthoPoly uses the default value of Nquad = 10 * (n + 1)
# I assign Nquad here so that I can print this value later
Nquad = 10 * (n + 1);
@time opq = OrthoPoly(
    "ex1_op",
    n,
    measure;
    addQuadrature = true,
)

# The PolyChaos function nw means nodes and weights.
# It returns the quadrature rule stored in opq in matrix form,
# with column 1 containing the nodes and column 2 containing the weights.
rule_polychaos = nw(opq)

row_count = size(rule_polychaos, 1)
nodes_polychaos = Vector{Float64}(undef, row_count)
weights_polychaos = Vector{Float64}(undef, row_count)
for i in 1:row_count
    nodes_polychaos[i] = rule_polychaos[i, 1]
    weights_polychaos[i] = rule_polychaos[i, 2]
end

@time nodes_momentdets_raw, weights_momentdets_raw =
    custom_gauss_quad_all_fn(moment_fn, which_f, n)
nodes_momentdets = Float64[Float64(x) for x in nodes_momentdets_raw]
weights_momentdets = Float64[Float64(w) for w in weights_momentdets_raw]

println("n = ", n)
println("Nquad = ", Nquad, "\n")

println("        nodes_polychaos                 nodes_momentdets")
for i in 1:n
    @printf("%2d   %.16e   %.16e\n", i, nodes_polychaos[i], nodes_momentdets[i])
end
print("\n")

println("       weights_polychaos               weights_momentdets")
for i in 1:n
    @printf("%2d   %.16e   %.16e\n", i, weights_polychaos[i], weights_momentdets[i])
end
print("\n")

println("maximum(abs.(nodes_polychaos .- nodes_momentdets)) = ",
    maximum(abs.(nodes_polychaos .- nodes_momentdets)))

println("maximum(abs.((weights_polychaos .- weights_momentdets) ./ weights_momentdets)) = ",
    maximum(abs.((weights_polychaos .- weights_momentdets) ./ weights_momentdets)), "\n")


# The Stieltjes comparison below is disabled in Example 1 because it takes too
# long for this weight. It is kept here so that it can be copied into a later
# example if needed.
#=
# Compare the same PolyChaos Gauss rule with the Gauss rule obtained from
# CustomGaussQuadrature using the Stieltjes procedure.

println("Compare the PolyChaos Gauss rule with")
println("CustomGaussQuadrature using the Stieltjes procedure.", "\n")

lnf_typed_fn = function (::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "exp(-abs(x)^k)"
    T_k = materialize_scalar_spec_fn(T, which_f[3])
    @assert T_k > zero(T)
    -(abs(T(x))^T_k)
end

mu0 = moment_fn(BigFloat, which_f, 0)
a, b = which_f[2]

@time begin
    a_vec_stieltjes, b_vec_stieltjes, nbits_stieltjes, r_stieltjes =
        stieltjes_a_vec_b_vec_final_fn(
            n,
            mu0,
            lnf_typed_fn,
            which_f,
            a,
            b,
        )

    nodes_stieltjes_raw, weights_stieltjes_raw =
        stieltjes_custom_gauss_quad_all_fn(
            n,
            mu0,
            a_vec_stieltjes,
            b_vec_stieltjes,
            a,
            b,
        )
    nodes_stieltjes = Float64[Float64(x) for x in nodes_stieltjes_raw]
    weights_stieltjes = Float64[Float64(w) for w in weights_stieltjes_raw]

    println("r_stieltjes = ", r_stieltjes)
    println("nbits_stieltjes = ", nbits_stieltjes, "\n")

    println("        nodes_polychaos      nodes_stieltjes")
    for i in 1:n
        @printf("%2d   %.16e   %.16e\n", i, nodes_polychaos[i], nodes_stieltjes[i])
    end
    print("\n")

    println("       weights_polychaos    weights_stieltjes_expabs")
    for i in 1:n
        @printf("%2d   %.16e   %.16e\n", i, weights_polychaos[i], weights_stieltjes[i])
    end
    print("\n")

    println("maximum(abs.(nodes_polychaos .- nodes_stieltjes)) = ",
        maximum(abs.(nodes_polychaos .- nodes_stieltjes)))

    println("maximum(abs.((weights_polychaos .- weights_stieltjes) ./ weights_stieltjes)) = ",
        maximum(abs.((weights_polychaos .- weights_stieltjes) ./ weights_stieltjes)), "\n")
end
=#
end


println("\n", "Example 2: weight is exp(-x³/3) for x > 0")
# Consider the weight function considered by Gautschi (1983),
# evaluated at x. This weight function is     
#   f(x) = exp(-x³/3) for x > 0
#        = 0  otherwise.
# Reference
# Gautschi, W. (1983) How and how not to check Gaussian quadrature
#formulae. BIT, 23, 209-216.      

# A let block gives this example its own local scope, as in Example 1.
let
which_f = ["chemistry example", [0, Inf]]
a, b = which_f[2]

n = 15

println("Compare the PolyChaos Gauss rule with")
println("CustomGaussQuadrature using moment determinants.", "\n")

# The chemistry-example moments are already stored in CustomGaussQuadrature:
#
#   μₛ = ∫₀∞ xˢ exp(-x³/3) dx
#      = 3^((s - 2)/3) Γ((s + 1)/3).
#
# Therefore we can use the stored moment function directly.
moment_fn = moment_stored_fn

measure = Measure("ex2_meas", x -> exp(-(x^3) / 3.0), (a, b), false)

Nquad = 10 * (n + 1);
@time opq = OrthoPoly(
    "ex2_op",
    n,
    measure;
    addQuadrature = true,
)

rule_polychaos = nw(opq)

row_count = size(rule_polychaos, 1)
nodes_polychaos = Vector{Float64}(undef, row_count)
weights_polychaos = Vector{Float64}(undef, row_count)
for i in 1:row_count
    nodes_polychaos[i] = rule_polychaos[i, 1]
    weights_polychaos[i] = rule_polychaos[i, 2]
end

@time nodes_momentdets_raw, weights_momentdets_raw =
    custom_gauss_quad_all_fn(moment_fn, which_f, n)
nodes_momentdets = Float64[Float64(x) for x in nodes_momentdets_raw]
weights_momentdets = Float64[Float64(w) for w in weights_momentdets_raw]

println("n = ", n)
println("Nquad = ", Nquad, "\n")

println("        nodes_polychaos                 nodes_momentdets")
for i in 1:n
    @printf("%2d   %.16e   %.16e\n", i, nodes_polychaos[i], nodes_momentdets[i])
end
print("\n")

println("       weights_polychaos               weights_momentdets")
for i in 1:n
    @printf("%2d   %.16e   %.16e\n", i, weights_polychaos[i], weights_momentdets[i])
end
print("\n")

println("maximum(abs.(nodes_polychaos .- nodes_momentdets)) = ",
    maximum(abs.(nodes_polychaos .- nodes_momentdets)))

println("maximum(abs.((weights_polychaos .- weights_momentdets) ./ weights_momentdets)) = ",
    maximum(abs.((weights_polychaos .- weights_momentdets) ./ weights_momentdets)), "\n")
end


println("\n", "Example 3: scaled chi pdf with m = 160")

# A let block gives this example its own local scope, as in Examples 1 and 2.
let
m = 160    
which_f = ["scaled chi pdf", [0,Inf], m]
a, b = which_f[2]

n = 33

println("Compare the PolyChaos Gauss rule with")
println("CustomGaussQuadrature using moment determinants.", "\n")

# The scaled-chi-pdf moments are already stored in CustomGaussQuadrature:
#
#   μₛ = ∫₀∞ xˢ f(x) dx
#      = (2/m)^(s/2) Γ((s + m)/2) / Γ(m/2).
#
# Therefore we can use the stored moment function directly.
moment_fn = moment_stored_fn

measure = Measure(
    "ex3_meas",
    x -> scaled_chi_pdf_fn(Float64, x, m),
    (a, b),
    false,
)

Nquad = 10 * (n + 1);
@time opq = OrthoPoly(
    "ex3_op",
    n,
    measure;
    addQuadrature = true,
)

rule_polychaos = nw(opq)

row_count = size(rule_polychaos, 1)
nodes_polychaos = Vector{Float64}(undef, row_count)
weights_polychaos = Vector{Float64}(undef, row_count)
for i in 1:row_count
    nodes_polychaos[i] = rule_polychaos[i, 1]
    weights_polychaos[i] = rule_polychaos[i, 2]
end

@time nodes_momentdets_raw, weights_momentdets_raw =
    custom_gauss_quad_all_fn(moment_fn, which_f, n)
nodes_momentdets = Float64[Float64(x) for x in nodes_momentdets_raw]
weights_momentdets = Float64[Float64(w) for w in weights_momentdets_raw]

println("n = ", n)
println("Nquad = ", Nquad, "\n")

println("        nodes_polychaos                 nodes_momentdets")
for i in 1:n
    @printf("%2d   %.16e   %.16e\n", i, nodes_polychaos[i], nodes_momentdets[i])
end
print("\n")

println("       weights_polychaos               weights_momentdets")
for i in 1:n
    @printf("%2d   %.16e   %.16e\n", i, weights_polychaos[i], weights_momentdets[i])
end
print("\n")

println("maximum(abs.(nodes_polychaos .- nodes_momentdets)) = ",
    maximum(abs.(nodes_polychaos .- nodes_momentdets)))

println("maximum(abs.((weights_polychaos .- weights_momentdets) ./ weights_momentdets)) = ",
    maximum(abs.((weights_polychaos .- weights_momentdets) ./ weights_momentdets)), "\n")
end
