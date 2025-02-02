module CustomGaussQuadrature
# Written by Dr. Paul Kabaila, 
# Department of Mathematical and Physical Sciences,
# La Trobe University, Melbourne, Australia.
# February 2025.
# 
# This module consists of the following Julia scripts:
#
# (1) gauss_quad_moment_dets_scr.jl
# This script is for computing the custom-made Gauss 
# quadrature nodes and weights using a moment-based method 
# via moment determinants to compute the coefficients
# in the three-term recurrence relation, using
# the type T, which is specified by the user. 
# This is followed by the use of the package
# GenericLinearAlgebra.jl
#
# (2) inner_prod_fns_scr.jl
# This script ...


using GaussQuadrature
using FastGaussQuadrature
using QuadGK
using SpecialFunctions
using OffsetArrays
using LinearAlgebra
using GenericLinearAlgebra
using DoubleFloats
using Distributions
using Polynomials
using Plots

# (1) 
include("gauss_quad_moment_dets_scr.jl")
export moment_fn
export μ_offsetvec_fn 
export Δ_fn
export Δ′_fn
export Δ_offsetvec_fn
export Δ′_offsetvec_fn
export α_offsetvec_fn
export β_vec_fn
export α_offsetvec_β_vec_fn
export a_vec_b_vec_μ₀_fn
export step1a_fn
export step1_fn
export step2_fn
export custom_gauss_quad_fn
export custom_gauss_quad_all_fn

# (2)
include("inner_prod_fns_scr.jl")
export ϕ_fn
export ϕ′_fn
export nodes_weights_support_fn
export scaled_chi_pdf_fn
export inner_prod_fns_fn
export ln_ϕ′_fn
export ln_scaled_chi_pdf_fn
export lnf_chemistry_fn
export lnf_hermite_fn
export lnf_laguerre_fn
export nodes_lnweights_support_fn
export inner_prod_fns_nonan_fn

# (3)
include("gauss_quad_stjieltjes_scr.jl")
export stjieltjes_a_vec_b_vec_nonan_fn
export stjieltjes_a_vec_b_vec_choosenbits_fn
export stjieltjes_a_vec_b_vec_final_fn
export stjieltjes_step2_fn
export stjieltjes_custom_gauss_quad_all_fn

# (4)
include("gauss_quad_utilities_scr.jl")
export integrand_transf_integral_fn
export plot_cdf_discrete_rv_fn
export scaled_chi_cdf_fn
export scaled_chi_cdf_fn


end
