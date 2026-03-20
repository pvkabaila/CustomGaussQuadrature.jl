module CustomGaussQuadrature
# Written by Dr. Paul Kabaila, 
# Department of Mathematical and Physical Sciences,
# La Trobe University, Melbourne, Australia.
# March 2026.
# 
# This module consists of the following Julia scripts:
#
# (1) gauss_quad_moment_dets_scr.jl
# This script is for computing the custom-made Gauss quadrature 
# rule with n nodes using a moment-based method via moment 
# determinants. This is done using the following two steps:
#
# Step 1: Compute the recursion coefficients 
# [α₀, α₁, ..., αₙ₋₁] and [β₁, ..., βₙ] of the 
# three-term recurrence relation using a moment-based 
# method via moment determinants. The aim of 
# Step 1 is to compute approximations to the 
# vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and 
# [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁] with maximum
# absolute errors and maximum relative errors,
# respectively, bounded above by 10⁻¹⁸.
# The map from the vector of moments
# [μ₀, μ₁, ..., μ₂ₙ₋₁] to the vector 
# [α₀, α₁, ..., αₙ₋₁,β₁, ..., βₙ] of recursion 
# coefficients is severely ill-conditioned. 
# This limitation is overcome by using BigFloat
# arithmetic with carefully chosen precision.
#
# Step 2: Compute the Gauss quadrature rule with n nodes
# from the eigenvalues and eigenvectors of the symmetric 
# tridiagonal Jacobi matrix Jₙ with main diagonal
# [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and subdiagonal 
# [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], using  
# Double64 arithmetic. The computation of these 
# eigenvalues and eigenvectors is carried out using
# the package GenericLinearAlgebra.jl
#
# (2) inner_prod_fns_scr.jl
# This script is for computing an approximation to the inner 
# product of two functions u: ℝ → ℝ and v: ℝ → ℝ. 
# Let f: ℝ → [0,∞) denote a specified nonnegative weight function. 
# The inner product of the two functions u and v is defined to be
#
#    ∞
#    ∫ u(x) v(x) f(x) dx.                                  
#   -∞
#
# To compute this inner product let the function g = uv, so that
# this integral becomes 
# 
#    ∞
#    ∫ g(x) f(x) dx.                                       
#   -∞
#
# We also suppose that the support of the weight function f
# is an interval with lower and upper endpoints a and b, 
# respectively. Here -∞ ≤ a < b ≤ ∞. Therefore this integral
# becomes
# 
#    b
#    ∫ g(x) f(x) dx.                                      (1) 
#    a
#
# To compute this integral, we use the transformation described
# on p.94 of Gautschi (2004). In other words, we transform the  
# support interval with lower and upper endpoints a and b, 
# respectively, to the interval [-1,1] using the 
# transformation
#
#   b                1
#   ∫ g(x) f(x) dx = ∫ g(ϕ(y)) f(ϕ(y)) ϕ′(y) dy              (2)
#   a               -1
#
# where  
#
#        ϕ(y) = (1/2) (b - a) y + (1/2) (b + a) if -∞ < a < b < ∞
#        
#        ϕ(y) = b - (1 - y) / (1 + y)           if -∞ = a < b < ∞
#   
#        ϕ(y) = a + (1 + y) / (1 - y)           if -∞ < a < b = ∞
#
#        ϕ(y) = y / (1 - y²)                    if -∞ = a < b = ∞
#
# The integral on the right-hand side of (1) is approximated using Gauss
# Legendre quadrature with r nodes. 
#
# (3) gauss_quad_stieltjes_scr.jl
# This script is for computing the custom-made Gauss quadrature 
# rule with n nodes using the Stieltjes procedure.
# This is done using the following two steps:
#
# Step 1: Compute the recursion coefficients 
# [α₀, α₁, ..., αₙ₋₁] and [β₁, ..., βₙ] of the 
# three-term recurrence relation using the 
# Stieltjes procedure, where the inner product 
# of two functions is found using the script 
# inner_prod_fns_scr.jl. The aim of 
# Step 1 is to compute approximations to the 
# vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and 
# [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁] with maximum
# absolute errors and maximum relative errors,
# respectively, bounded above by 10⁻¹⁸.
#
# Step 2: Compute the Gauss quadrature rule with n nodes
# from the eigenvalues and eigenvectors of the symmetric 
# tridiagonal Jacobi matrix Jₙ with main diagonal
# [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and subdiagonal 
# [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], using  
# Double64 arithmetic. The computation of these 
# eigenvalues and eigenvectors is carried out using
# the package GenericLinearAlgebra.jl
#
# Reference
# Gautschi, W. (2004) Orthogonal Polynomials, Computation and 
# Approximation. Oxford University Press, Oxford.

using FastGaussQuadrature
using GaussQuadrature
using SpecialFunctions
using OffsetArrays
using LinearAlgebra
using GenericLinearAlgebra
using DoubleFloats
using Polynomials

# (1) 
include("gauss_quad_moment_dets_scr.jl")
export moment_stored_fn
export custom_gauss_quad_fn
export custom_gauss_quad_all_fn

# (2)
include("inner_prod_fns_scr.jl")
export scaled_chi_pdf_fn
export ln_scaled_chi_pdf_fn
export lnf_chemistry_fn
export lnf_hermite_fn
export lnf_laguerre_fn

# (3)
include("gauss_quad_stieltjes_scr.jl")
export stieltjes_a_vec_b_vec_nonan_fn
export stieltjes_a_vec_b_vec_choosenbits_fn
export stieltjes_a_vec_b_vec_final_fn
export stieltjes_step2_fn
export stieltjes_custom_gauss_quad_all_fn

end
