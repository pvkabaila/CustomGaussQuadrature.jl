# R interface script for Method 2 (Stieltjes procedure)
# with a built-in (stored) weight function.
#
# PURPOSE
# -------
# This script is called from R via JuliaConnectoR to compute
# custom Gauss quadrature nodes and weights using the Stieltjes
# procedure with one of the built-in weight functions
# (those for which the log-weight function is already stored in the package).
#
# PRE-CONDITIONS
# --------------
# Before including this script the following must already be defined
# in the Julia session:
#
#   using CustomGaussQuadrature
#
#   which_f -- Array identifying a built-in weight function, e.g.
#              ["scaled chi pdf", [0, Inf], 160]
#   n       -- Number of Gauss quadrature nodes (positive integer).
#
# No user-supplied .jl file is needed for the weight function itself.
# which_f and n can be defined directly inside a juliaEval() call.
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
#   juliaEval('which_f = ["scaled chi pdf", [0, Inf], 160]; n = 33')
#   juliaEval('include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))),
#       "src", "r_interface_stieltjes_stored_scr.jl"))')
#   nodes   <- juliaEval('cgq_nodes')
#   weights <- juliaEval('cgq_weights')

include(joinpath(@__DIR__, "stieltjes_lnf_stored_scr.jl"))
cgq_nodes   = convert(Vector{Float64}, stieltjes_nodes)
cgq_weights = convert(Vector{Float64}, stieltjes_weights)
