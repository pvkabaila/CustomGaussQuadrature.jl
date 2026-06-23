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
#   A Julia function with signature
#       weight_name_moment_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
#   returning the s'th moment of the weight function.
#   Here weight_name_moment_fn is derived from which_f[1] by replacing
#   each run of non-alphanumeric characters by "_" and then appending
#   "_moment_fn". For example,
#       "weibull pdf"    -> weibull_pdf_moment_fn
#   which_f    -- Array identifying the weight function, e.g.
#                 ["weibull pdf", [0, Inf], "2.1"]
#   n          -- Number of Gauss quadrature nodes (positive integer).
#
# These must be defined in the Julia session before this script is
# included. In the current R workflow this naming convention is handled
# by an R helper routine, so the R user only needs to provide which_f,
# n and the body of the moment function; see the typical R usage below.
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
#   juliaEval('using Pkg')
#   juliaEval('Pkg.develop(path=raw"C:/path/to/CustomGaussQuadrature")')
#   juliaEval('using CustomGaussQuadrature')
#   juliaEval('pathof(CustomGaussQuadrature)')
#   juliaEval('using SpecialFunctions')
#   # Use strings for finite non-integer constants so Julia preserves
#   # the decimal specification. Vector or list parameters may contain
#   # these scalar specifications recursively.
#   # Example 1: scalar parameter specification.
#   which_f <- list("weibull pdf", c(0, Inf), "2.1")
#   moment_fn_body <- '
#       @assert which_f[1] == "weibull pdf"
#       T_k = materialize_scalar_spec_fn(T, which_f[3])
#       @assert T_k > 0
#       ...
#   '
#
#   # Example 2: vector-valued parameter specification.
#   which_f <- list("inverse gamma pdf", c(0, Inf), list("18.5", "3.2"))
#   moment_fn_body <- '
#       @assert which_f[1] == "inverse gamma pdf"
#       alpha_spec, beta_spec = which_f[3]
#       T_alpha = materialize_scalar_spec_fn(T, alpha_spec)
#       T_beta = materialize_scalar_spec_fn(T, beta_spec)
#       @assert T_alpha > convert(T, s)
#       ...
#   '
#   n <- 9
#   run_moment_dets_example <- function(which_f, n, moment_fn_body) {
#     julia_value_literal <- function(value) {
#       if (is.character(value) && length(value) == 1) {
#         escaped <- gsub("\\\\", "\\\\\\\\", value)
#         escaped <- gsub('"', '\\"', escaped, fixed = TRUE)
#         return(sprintf('"%s"', escaped))
#       }
#       if (is.numeric(value) && length(value) == 1) {
#         if (is.infinite(value)) {
#           return(if (value > 0) "Inf" else "-Inf")
#         }
#         return(format(value, trim = TRUE, scientific = FALSE))
#       }
#       if (is.list(value) || (is.atomic(value) && length(value) > 1)) {
#         item_literals <- vapply(value, julia_value_literal, character(1))
#         return(sprintf("[%s]", paste(item_literals, collapse = ", ")))
#       }
#       stop("This example supports scalar numeric or character specifications, or vectors/lists containing them.")
#     }
#
#     weight_name <- which_f[[1]]
#     support <- which_f[[2]]
#     left_endpoint <- if (is.infinite(support[1])) "-Inf" else format(support[1], trim = TRUE)
#     right_endpoint <- if (is.infinite(support[2])) "Inf" else format(support[2], trim = TRUE)
#     support_str <- paste0(left_endpoint, ", ", right_endpoint)
#
#     if (length(which_f) == 2) {
#       juliaEval(sprintf('which_f = ["%s", [%s]]', weight_name, support_str))
#     } else {
#       parameter_str <- julia_value_literal(which_f[[3]])
#       juliaEval(sprintf('which_f = ["%s", [%s], %s]', weight_name, support_str, parameter_str))
#     }
#     juliaEval(sprintf('n = %d', n))
#
#     fn_name <- paste0(gsub("[^A-Za-z0-9]+", "_", weight_name), "_moment_fn")
#     juliaEval(sprintf(
# 'function %s(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
#     @assert s >= 0
# %s
# end',
#       fn_name, moment_fn_body))
#
#     juliaEval('include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))), "src", "r_interface_moment_dets_new_scr.jl"))')
#
#     list(
#       nodes = juliaEval('cgq_nodes'),
#       weights = juliaEval('cgq_weights')
#     )
#   }
#
#   result <- run_moment_dets_example(which_f, n, moment_fn_body)
#   nodes   <- result$nodes
#   weights <- result$weights
#

moment_fn_sym = Symbol(replace(which_f[1], r"[^A-Za-z0-9]+" => "_") * "_moment_fn")
@assert isdefined(Main, moment_fn_sym)
moment_fn = getfield(Main, moment_fn_sym)
cgq_nodes_raw, cgq_weights_raw = custom_gauss_quad_all_fn(moment_fn, which_f, n)
cgq_nodes   = convert(Vector{Float64}, cgq_nodes_raw)
cgq_weights = convert(Vector{Float64}, cgq_weights_raw)
