# R interface script for Method 2 (Stieltjes procedure)
# with a user-defined log-weight function.
#
# PURPOSE
# -------
# This script is called from R via JuliaConnectoR to compute
# custom Gauss quadrature nodes and weights using the Stieltjes
# procedure with a user-defined log-weight function.
#
# PRE-CONDITIONS
# --------------
# Before including this script the following must already be defined
# in the Julia session:
#
#   using CustomGaussQuadrature
#
#   T       -- floating-point type, e.g. T = BigFloat
#   which_f -- Array identifying the weight function, e.g.
#              ["weibull pdf", [0, Inf], 2.0]
#   n       -- Number of Gauss quadrature nodes (positive integer).
#   lnf_fn  -- Julia function of x evaluating log(f(x)), e.g.
#              lnf_fn = x -> lnf_weibull_pdf_fn(T, which_f, x)
#   mu0     -- Zeroth moment: integral of f over the support.
#              For a pdf, mu0 = convert(T, 1).
#
# These must be defined by including a user-supplied .jl file, e.g.:
#   juliaEval('include("C:/path/to/my_lnf_fn.jl")')
#
# POST-CONDITIONS
# ---------------
# After this script completes, the following are defined:
#
#   cgq_nodes   -- Vector{Float64} of n Gauss quadrature nodes.
#   cgq_weights -- Vector{Float64} of n Gauss quadrature weights.
#
# TYPICAL R USAGE
# ---------------
#   library(JuliaConnectoR)
#   juliaEval('using CustomGaussQuadrature')
#   juliaEval('include("C:/path/to/my_lnf_fn.jl")')
#   juliaEval('include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))),
#       "src", "r_interface_stieltjes_new_scr.jl"))')
#   nodes   <- juliaEval('cgq_nodes')
#   weights <- juliaEval('cgq_weights')

include(joinpath(@__DIR__, "stieltjes_lnf_new_scr.jl"))
cgq_nodes   = convert(Vector{Float64}, stieltjes_nodes)
cgq_weights = convert(Vector{Float64}, stieltjes_weights)
