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
#              ["weibull pdf", [0, Inf], "2.1"]
#   n       -- Number of Gauss quadrature nodes (positive integer).
#   lnf_typed_fn -- Julia function of (T, which_f, x) evaluating log(f(x)), e.g.
#                   lnf_typed_fn = lnf_weibull_pdf_fn
#   mu0     -- Zeroth moment: integral of f over the support.
#              For a pdf, use mu0 = 1. Use a string only for a finite
#              non-integer constant whose decimal value matters.
#   j_max   -- Optional positive integer overriding the default cap used
#              when increasing r in the Stieltjes procedure.
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
#   cgq_r       -- Integer giving the final support-quadrature node count r
#                  chosen by the Stieltjes driver.
#
# TYPICAL R USAGE
# ---------------
#   library(JuliaConnectoR)
#   juliaEval('Pkg.develop(path=raw"C:/path/to/CustomGaussQuadrature")')
#   juliaEval('using CustomGaussQuadrature')
#   juliaEval('pathof(CustomGaussQuadrature)')
#   mu0_body <- '1'
#   # Use strings for finite non-integer constants so Julia preserves
#   # the decimal specification. Vector or list parameters may contain
#   # these scalar specifications recursively. mu0_body remains a Julia
#   # source fragment.
#   # Example 1: scalar parameter specification.
#   which_f <- list("weibull pdf", c(0, Inf), "2.1")
#   lnf_typed_fn_body <- '
#       @assert which_f[1] == "weibull pdf"
#       T_k = materialize_scalar_spec_fn(T, which_f[3])
#       @assert T_k > zero(T_k)
#       ...
#   '
#
#   # Example 2: vector-valued parameter specification.
#   which_f <- list("inverse gamma pdf", c(0, Inf), list("18.5", "3.2"))
#   lnf_typed_fn_body <- '
#       @assert which_f[1] == "inverse gamma pdf"
#       alpha_spec, beta_spec = which_f[3]
#       T_alpha = materialize_scalar_spec_fn(T, alpha_spec)
#       T_beta = materialize_scalar_spec_fn(T, beta_spec)
#       ...
#   '
#   n <- 9
#   julia_value_literal <- function(value) {
#     if (is.character(value) && length(value) == 1) {
#       escaped <- gsub("\\\\", "\\\\\\\\", value)
#       escaped <- gsub('"', '\\"', escaped, fixed = TRUE)
#       return(sprintf('"%s"', escaped))
#     }
#     if (is.numeric(value) && length(value) == 1) {
#       if (is.infinite(value)) {
#         return(if (value > 0) "Inf" else "-Inf")
#       }
#       return(format(value, trim = TRUE, scientific = FALSE))
#     }
#     if (is.list(value) || (is.atomic(value) && length(value) > 1)) {
#       item_literals <- vapply(value, julia_value_literal, character(1))
#       return(sprintf("[%s]", paste(item_literals, collapse = ", ")))
#     }
#     stop("This example supports scalar numeric or character specifications, or vectors/lists containing them.")
#   }
#
#   run_stieltjes_new_example <- function(which_f, n, lnf_typed_fn_body, mu0_body, j_max = 40) {
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
#     fn_name <- paste0(gsub("[^A-Za-z0-9]+", "_", weight_name), "_lnf_fn")
#     juliaEval(sprintf(
# 'function %s(T::Type, which_f, x::AbstractFloat)
# %s
# end',
#       fn_name, lnf_typed_fn_body))
#     # The R helper defines only the three-argument function
#     # lnf_typed_fn(T, which_f, x). It does NOT build the final
#     # one-argument closure lnf_fn(x). That closure is built later,
#     # inside the Stieltjes driver, after the driver has chosen the
#     # working type T for the current branch.
#     juliaEval(sprintf('lnf_typed_fn = %s', fn_name))
#     juliaEval(sprintf('mu0 = %s', mu0_body))
#     juliaEval(sprintf('j_max = %d', j_max))
#
#     juliaEval('include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))),
#         "src", "r_interface_stieltjes_new_scr.jl"))')
#
#     list(
#       nodes = juliaEval('cgq_nodes'),
#       weights = juliaEval('cgq_weights'),
#       r = juliaEval('cgq_r')
#     )
#   }
#
#   result <- run_stieltjes_new_example(which_f, n, lnf_typed_fn_body, mu0_body, j_max = 80)
#   nodes   <- result$nodes
#   weights <- result$weights
#   r       <- result$r
#
#   The next helper profiles the main R-to-Julia stages separately.
#   It times each juliaEval call used to set up the Julia variables,
#   define lnf_typed_fn, run the include that performs the Stieltjes
#   computation, and fetch the nodes and weights back into R.
#
#   This is more informative than applying system.time() only to the
#   whole call run_stieltjes_new_example(...). In that overall timing,
#   the reported R user and system CPU times are often very small because
#   R is mostly waiting while JuliaConnectoR hands work to Julia and while
#   Julia performs startup, JIT compilation, and the numerical work.
#   Consequently the elapsed time can be several seconds even though the
#   R CPU times are near zero, so a single system.time(...) does not show
#   which stage is actually responsible for the delay.
#
#   profile_run_stieltjes_new_example <- function(which_f, n, lnf_typed_fn_body, mu0_body, j_max = 40) {
#     timings <- data.frame(
#       step = character(),
#       user = numeric(),
#       system = numeric(),
#       elapsed = numeric(),
#       stringsAsFactors = FALSE
#     )
#
#     capture_step <- function(label, expr) {
#       tm <- system.time(value <- eval(substitute(expr), envir = parent.frame()))
#       timings <<- rbind(
#         timings,
#         data.frame(
#           step = label,
#           user = unname(tm[["user.self"]]),
#           system = unname(tm[["sys.self"]]),
#           elapsed = unname(tm[["elapsed"]]),
#           stringsAsFactors = FALSE
#         )
#       )
#       value
#     }
#
#     weight_name <- which_f[[1]]
#     support <- which_f[[2]]
#     left_endpoint <- if (is.infinite(support[1])) "-Inf" else format(support[1], trim = TRUE)
#     right_endpoint <- if (is.infinite(support[2])) "Inf" else format(support[2], trim = TRUE)
#     support_str <- paste0(left_endpoint, ", ", right_endpoint)
#
#     if (length(which_f) == 2) {
#       capture_step("set which_f",
#         juliaEval(sprintf('which_f = ["%s", [%s]]', weight_name, support_str))
#       )
#     } else {
#       parameter_str <- julia_value_literal(which_f[[3]])
#       capture_step("set which_f",
#         juliaEval(sprintf('which_f = ["%s", [%s], %s]', weight_name, support_str, parameter_str))
#       )
#     }
#
#     capture_step("set n",
#       juliaEval(sprintf('n = %d', n))
#     )
#
#     fn_name <- paste0(gsub("[^A-Za-z0-9]+", "_", weight_name), "_lnf_fn")
#     capture_step("define lnf fn",
#       juliaEval(sprintf(
# 'function %s(T::Type, which_f, x::AbstractFloat)
# %s
# end',
#         fn_name, lnf_typed_fn_body))
#     )
#
#     capture_step("bind lnf fn",
#       juliaEval(sprintf('lnf_typed_fn = %s', fn_name))
#     )
#
#     capture_step("set mu0",
#       juliaEval(sprintf('mu0 = %s', mu0_body))
#     )
#
#     capture_step("set j_max",
#       juliaEval(sprintf('j_max = %d', j_max))
#     )
#
#     capture_step("run include",
#       juliaEval('t_include = @elapsed include(joinpath(dirname(dirname(pathof(CustomGaussQuadrature))), "src", "r_interface_stieltjes_new_scr.jl"))')
#     )
#
#     julia_include_elapsed <- capture_step("get Julia elapsed",
#       juliaEval('t_include')
#     )
#
#     nodes <- capture_step("get nodes",
#       juliaEval('cgq_nodes')
#     )
#
#     weights <- capture_step("get weights",
#       juliaEval('cgq_weights')
#     )
#
#     r <- capture_step("get r",
#       juliaEval('cgq_r')
#     )
#
#     print(timings, row.names = FALSE)
#     cat(sprintf("Julia include elapsed: %.4f seconds\n", julia_include_elapsed))
#
#     list(
#       result = list(nodes = nodes, weights = weights, r = r),
#       timings = timings,
#       julia_include_elapsed = julia_include_elapsed
#     )
#   }
#
#   profile <- profile_run_stieltjes_new_example(which_f, n, lnf_typed_fn_body, mu0_body, j_max = 80)
#   nodes   <- profile$result$nodes
#   weights <- profile$result$weights
#   r       <- profile$result$r

include(joinpath(@__DIR__, "stieltjes_lnf_new_scr.jl"))
cgq_nodes   = convert(Vector{Float64}, stieltjes_nodes)
cgq_weights = convert(Vector{Float64}, stieltjes_weights)
cgq_r       = r
