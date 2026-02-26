# The latest version of my Julia package CustomGaussQuadrature is in the folder:
# \RESEARCH - NUMERICAL METHODS\QUADRATURE\Custom GAUSS\CustomGaussQuadrature - Julia package\
# CustomGaussQuadrature\

# This script tests the version on my computer, not the
# version in the Julia General Registry.

# To run the code in this package, the first step in VS code is 
# File > Open Folder... > open the above folder

# Use CTRL + Shift + P to get to the Command Palette in VS code and choose
# Julia: Start REPL
# to start the Julia REPL.

# Use 
# Go to package mode by typing ]
# Activate the package environment by running the 
# activate . command
# to activate the package. Exit package mode using Ctrl + C

# Then use the following command to run this script:
# julia> include("test/test_stjieltjes_lnf_new_scr.jl")

# Copy the output at the REPL into a 
# text document. Do NOT execute this 
# code line-by-line in the REPL, as this introduces
# ugly extra spaces.

# Have a look at the README.md file by opening this file in VS code and then using
# CTRL + K      V 
# to get a preview of this file.

println("Test of the version of the CustomGaussQuadrature package")
println("on my computer, not the version in the Julia General Registry.")

using CustomGaussQuadrature
# This script tests the version on my computer, not the
# version in the Julia General Registry. See path.
pathof(CustomGaussQuadrature)

using Printf
using Plots
using GaussQuadrature
using SpecialFunctions

println("\n", "------------------------------------------------------")
#---------------------------------------
# Consider the "new" weight function 
# which is the Weibull pdf with scale 
# parameter λ=1.
#
# The Weibull pdf is 
# f(x; λ, k) = (k/λ) (x/λ)ᵏ⁻¹ exp(-(x/λ)ᵏ) for x ≥ 0
#            = 0                           otherwise
# where k > 0 is the shape parameter and λ > 0 is the
# scale parameter. The r'th raw moment of a random 
# variable X with this pdf is 
#      E(Xʳ) = λʳ Γ(1 + (r/k))
#
# We consider the particular case that the scale 
# parameter λ = 1. In this case the pdf is
# f(x; k) = k xᵏ⁻¹ exp(-xᵏ) for x ≥ 0
#         = 0               otherwise
# where k > 0 is the shape parameter.
# The r'th raw moment of a random 
# variable X with this pdf is 
#      E(Xʳ) = Γ(1 + (r/k))

using SpecialFunctions

function lnf_weibull_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    @assert x > convert(T,0)
    k = which_f[3]
    T_k = parse(T, string(k))
    @assert T_k > convert(T, 0)
    log(T_k) + (T_k - convert(T,1)) * log(x) - x^T_k
end

function μ₀_μ₁_weibull_pdf_fn(::Type{T}, which_f) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    k = which_f[3]
    T_k = parse(T, string(k))
    @assert T_k > convert(T, 0)
    T_1 = convert(T, 1)
    μ₀ = convert(T,1)
    μ₁ = gamma(T_1 + (T_1/T_k))
    [μ₀ μ₁]
end


k = 2.0
which_f = ["weibull pdf", [0, Inf], k]
n = 10

println("which_f = ", which_f)
println("n = ", n, "\n")

T = BigFloat

lnf_fn = x -> lnf_weibull_pdf_fn(T, which_f, x);
μ₀, μ₁ = μ₀_μ₁_weibull_pdf_fn(T, which_f);


#  @__DIR__ gives the folder of the file currently being executed
#  .. = parent folder (package root),
#  then into src,
#  then to stjieltjes_lnf_stored_scr.jl

#  The following command results in values for 
#  stjieltjes_nodes and stjieltjes_weights
include(joinpath(@__DIR__, "..", "src", "stjieltjes_lnf_new_scr.jl"));

using SpecialFunctions

function moment_weibull_pdf_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    k = which_f[3]
    @assert k > 0
    T_k = parse(T, string(k))
    @assert r ≥ 0
    if r == 0
      return(convert(T, 1))
    end
    T_1 = convert(T, 1)
    T_r = convert(T, r)
    gamma(T_1 + (T_r/T_k))
end

#-----------------------------------------------------------------
# Compute the custom Gauss quadrature nodes and weights
# using the moment determinant method
moment_fn = moment_weibull_pdf_fn
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

# Compare the nodes and weights computed by the Stjieltjes
# procedure with the nodes and weights computed using the
# moment determinant method
diff_nodes = stjieltjes_nodes - nodes;
rel_diff_weights = (stjieltjes_weights - weights) ./ weights;

diff_nodes = convert(Vector{Float64}, diff_nodes);
rel_diff_weights = convert(Vector{Float64}, rel_diff_weights);

println("      stjieltjes_nodes - nodes    (stjieltjes_weights - weights)./weights")
for i in 1:lastindex(diff_nodes)
	@printf "%2d     " i
	@printf "%.16e     " diff_nodes[i]
    @printf "%.16e  \n" rel_diff_weights[i]
end