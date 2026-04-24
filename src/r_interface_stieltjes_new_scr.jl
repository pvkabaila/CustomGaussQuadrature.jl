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
#   which_f -- Array identifying the weight function, e.g.
#              ["weibull pdf", [0, Inf], 2.0]
#   n       -- Number of Gauss quadrature nodes (positive integer).
#   lnf_user_fn -- Julia function of (T, which_f, x) evaluating log(f(x)), e.g.
#                  lnf_user_fn = lnf_weibull_pdf_fn
#   mu0     -- Zeroth moment: integral of f over the support.
#              For a pdf, mu0 = convert(Double64, 1).
#
# These must be defined in the Julia session before this script is
# included. In the current R workflow this setup is handled by an R
# helper routine, so the R user only needs to provide which_f, n, the
# body of the log-weight function and the expression for mu0; see the
# typical R usage below.
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
#   which_f <- list("weibull pdf", c(0, Inf), 2.1)
#   lnf_fn_body <- '
#       @assert which_f[1] == "weibull pdf"
#       ...
#   '
#   mu0_body <- 'convert(Double64, 1)'
#   n <- 9
#   run_stieltjes_new_example <- function(which_f, n, lnf_fn_body, mu0_body) {
#     weight_name <- which_f[[1]]
#     support <- which_f[[2]]
#     left_endpoint <- if (is.infinite(support[1])) "-Inf" else format(support[1], trim = TRUE)
#     right_endpoint <- if (is.infinite(support[2])) "Inf" else format(support[2], trim = TRUE)
#     support_str <- paste0(left_endpoint, ", ", right_endpoint)
#
#     if (length(which_f) == 2) {
#       juliaEval(sprintf('which_f = ["%s", [%s]]', weight_name, support_str))
#     } else {
#       juliaEval(sprintf('which_f = ["%s", [%s], %s]', weight_name, support_str, which_f[[3]]))
#     }
#     juliaEval(sprintf('n = %d', n))
#
#     fn_name <- paste0(gsub("[^A-Za-z0-9]+", "_", weight_name), "_lnf_fn")
#     juliaEval(sprintf(
# 'function %s(T::Type, which_f, x::AbstractFloat)
# %s
# end',
#       fn_name, lnf_fn_body))
#     juliaEval(sprintf('lnf_user_fn = %s', fn_name))
#     juliaEval(sprintf('mu0 = %s', mu0_body))
#
#     juliaEval('include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))),
#         "src", "r_interface_stieltjes_new_scr.jl"))')
#
#     list(
#       nodes = juliaEval('cgq_nodes'),
#       weights = juliaEval('cgq_weights')
#     )
#   }
#
#   result <- run_stieltjes_new_example(which_f, n, lnf_fn_body, mu0_body)
#   nodes   <- result$nodes
#   weights <- result$weights

include(joinpath(@__DIR__, "stieltjes_lnf_new_scr.jl"))
cgq_nodes   = convert(Vector{Float64}, stieltjes_nodes)
cgq_weights = convert(Vector{Float64}, stieltjes_weights)
