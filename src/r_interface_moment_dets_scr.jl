# R interface script for Method 1 (moment determinants)
# with a user-defined moment function.
#
# PURPOSE
# -------
# This script is called from R via JuliaConnectoR to compute
# custom Gauss quadrature nodes and weights using the moment
# determinants method with a user-defined moment function.
#
# PRE-CONDITIONS
# --------------
# Before including this script the following must already be defined
# in the Julia session:
#
#   using CustomGaussQuadrature
#
#   moment_fn  -- Julia function with signature
#                 moment_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
#                 returning the s'th moment of the weight function.
#   which_f    -- Array identifying the weight function, e.g.
#                 ["weibull pdf", [0, Inf], 2.0]
#   n          -- Number of Gauss quadrature nodes (positive integer).
#
# These must be defined by including a user-supplied .jl file, e.g.:
#   juliaEval('include("C:/path/to/my_moment_fn.jl")')
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
#   juliaEval('include("C:/path/to/my_moment_fn.jl")')
#   juliaEval('include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))),
#       "src", "r_interface_moment_dets_scr.jl"))')
#   nodes   <- juliaEval('cgq_nodes')
#   weights <- juliaEval('cgq_weights')

cgq_nodes_raw, cgq_weights_raw = custom_gauss_quad_all_fn(moment_fn, which_f, n)
cgq_nodes   = convert(Vector{Float64}, cgq_nodes_raw)
cgq_weights = convert(Vector{Float64}, cgq_weights_raw)
