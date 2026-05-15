# choice_of_r_1_scr.jl
# This script is for investigation, not for automated testing.

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
using Printf
using DoubleFloats

println("For the specified weight function, inspect how the Stieltjes driver")
println("increases r and how quickly the recurrence coefficients converge.", "\n")

println("\n", "Use plain reassignment of offset, j_max and epsilon before each experiment block.")

function inspect_r_search_fn(n, μ₀, lnf_typed_fn, which_f; offset=7, j_max=40, epsilon=1.0e-15)
    println("which_f = ", which_f)
    println("n = ", n, ", j_max = ", j_max, ", epsilon = ", epsilon)

    a, b = which_f[2]

    try
        j = 3
        r = j * (offset + n)
        stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits =
            CustomGaussQuadrature.stieltjes_a_vec_b_vec_choosenbits_fn(
                n,
                μ₀,
                lnf_typed_fn,
                which_f,
                a,
                b,
                r,
            )
        j = j + 1
        r = j * (offset + n)
        stieltjes_a_vec_new, stieltjes_b_vec_new, stieltjes_nbits_new =
            CustomGaussQuadrature.stieltjes_a_vec_b_vec_choosenbits_fn(
                n,
                μ₀,
                lnf_typed_fn,
                which_f,
                a,
                b,
                r,
            )

        abs_error_a = convert(Vector{Float64}, abs.(stieltjes_a_vec_new - stieltjes_a_vec))
        max_abs_error_a = maximum(abs_error_a)
        rel_error_b = convert(
            Vector{Float64},
            abs.((stieltjes_b_vec_new - stieltjes_b_vec) ./ stieltjes_b_vec_new),
        )
        max_rel_error_b = maximum(rel_error_b)

        println(
            "r search: j = ", j,
            ", r = ", r,
            ", nbits = ", stieltjes_nbits_new,
            ", max_abs_error_a = ", @sprintf("%.3e", max_abs_error_a),
            ", max_rel_error_b = ", @sprintf("%.3e", max_rel_error_b),
        )

        while (max_abs_error_a > epsilon || max_rel_error_b > epsilon)
            stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits =
                stieltjes_a_vec_new, stieltjes_b_vec_new, stieltjes_nbits_new
            if j == j_max
                throw(DomainError(j, " needed value of the r search index exceeds j_max"))
            end
            j = j + 1
            r = j * (offset + n)
            stieltjes_a_vec_new, stieltjes_b_vec_new, stieltjes_nbits_new =
                CustomGaussQuadrature.stieltjes_a_vec_b_vec_choosenbits_fn(
                    n,
                    μ₀,
                    lnf_typed_fn,
                    which_f,
                    a,
                    b,
                    r,
                )

            abs_error_a = convert(Vector{Float64}, abs.(stieltjes_a_vec_new - stieltjes_a_vec))
            max_abs_error_a = maximum(abs_error_a)
            rel_error_b = convert(
                Vector{Float64},
                abs.((stieltjes_b_vec_new - stieltjes_b_vec) ./ stieltjes_b_vec_new),
            )
            max_rel_error_b = maximum(rel_error_b)

            println(
                "r search: j = ", j,
                ", r = ", r,
                ", nbits = ", stieltjes_nbits_new,
                ", max_abs_error_a = ", @sprintf("%.3e", max_abs_error_a),
                ", max_rel_error_b = ", @sprintf("%.3e", max_rel_error_b),
            )
        end

        setprecision(BigFloat, 256, base=2)
        stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits =
            stieltjes_a_vec_new, stieltjes_b_vec_new, stieltjes_nbits_new
        println("r = ", r, ", nbits = ", stieltjes_nbits)
    catch err
        println("r search failed with ", typeof(err), ": ", err)
    end
    println()
end


function lnf_weibull_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    @assert x > convert(T,0)
    T_k = materialize_scalar_spec_fn(T, which_f[3])
    @assert T_k > convert(T, 0)
    log(T_k) + (T_k - convert(T,1)) * log(x) - x^T_k
end


println("Inspect the Stieltjes r search for a Weibull pdf weight function.")
    offset = 7;
    j_max = 40;
    epsilon = 1.0e-15;
    println("offset = ", offset, ", j_max = ", j_max, ", epsilon = ", epsilon)
     mu0 = 1
lnf_typed_fn = lnf_weibull_pdf_fn
for n in [6]    
    println("n = ", n)
    for k in ["3.1"]
        inspect_r_search_fn(
            n,
            mu0,
            lnf_typed_fn,
            ["weibull pdf", [0, Inf], k];
            offset=offset,
            j_max=j_max,
            epsilon=epsilon,
        )
    end
end

println("\n", "------------------------------------------------------")
println("Inspect the Stieltjes r search for a scaled chi pdf weight function.")
    offset = 7
    j_max = 40
    epsilon = 1.0e-15
    println("offset = ", offset, ", j_max = ", j_max, ", epsilon = ", epsilon)
    n = 10
    println("n = ", n)
for m in [1, 2, 5, 20, 160]
inspect_r_search_fn(
    n,
    convert(T, 1),
    (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3]),
    ["scaled chi pdf", [0, Inf], m];
        offset=offset,
        j_max=j_max,
    epsilon=epsilon,
)
end


println("\n", "------------------------------------------------------")
println("\n", "------------------------------------------------------")
println("For k = 1, the Weibull pdf weight function that we consider")
println("is the same as the Laguerre weight function as defined in") 
println("the Julia package GaussQuadrature.jl, with parameter α = 0.")

k = 1.0
which_f = ["weibull pdf", [0, Inf], k]
n = 10

println("which_f = ", which_f)
println("n = ", n, "\n")

lnf_typed_fn = lnf_weibull_pdf_fn;
mu0 = 1;

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

println("r =", r, "\n")

T = BigFloat

nodes, weights = laguerre(n, convert(T,0))

diff_nodes = stieltjes_nodes - nodes;
rel_diff_weights = (stieltjes_weights - weights) ./ weights;

diff_nodes = convert(Vector{Float64}, diff_nodes);
rel_diff_weights = convert(Vector{Float64}, rel_diff_weights);

println("      stieltjes_nodes - nodes    (stieltjes_weights - weights)./weights")
for i in 1:lastindex(diff_nodes)
	@printf "%2d     " i
	@printf "%.16e     " diff_nodes[i]
    @printf "%.16e  \n" rel_diff_weights[i]
end

print("\n")
println("maximum(abs.(stieltjes_nodes - nodes)) = ", 
convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((stieltjes_weights - weights) ./ weights)) = ", 
convert(Float64, maximum(abs.(rel_diff_weights))))