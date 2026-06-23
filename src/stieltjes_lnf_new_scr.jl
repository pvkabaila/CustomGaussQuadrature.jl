# User-defined Stieltjes driver.
#
# This file is intentionally not included from
# CustomGaussQuadrature.jl. It is meant to be run explicitly with
# include(...), typically from the REPL, tests, or an interface script.
#
# Required inputs in the caller:
#   which_f
#   n
#   lnf_typed_fn
#   mu0
# Optional inputs:
#   offset, j_max, epsilon
#
# Outputs in the caller:
#   nodes_stieltjes
#   weights_stieltjes
#   a_vec_stieltjes
#   b_vec_stieltjes
#   nbits_stieltjes
#   r
#
# Backward-compatible aliases also provided in the caller:
#   stieltjes_nodes
#   stieltjes_weights
#   stieltjes_a_vec
#   stieltjes_b_vec
#   stieltjes_nbits

using CustomGaussQuadrature: Double64

a, b = which_f[2];
# Keep the support endpoints in their input representation.
# The Stieltjes driver chooses the working type T internally, and the
# support-quadrature helpers convert a and b to that local T when needed.
# This prevents the stored type of the endpoints from forcing BigFloat or
# Double64 arithmetic before the driver has chosen which branch to test.

if @isdefined mu0
    @assert mu0 isa AbstractString || mu0 isa Integer
    μ₀_input = mu0
else
    throw(UndefVarError(:mu0))
end

μ₀_double64 = CustomGaussQuadrature.stieltjes_materialize_typed_scalar_fn(
    μ₀_input,
    Double64,
    which_f,
)

# This script keeps lnf_typed_fn and which_f separate and passes both to the
# driver. The driver then chooses T and creates the final one-argument lnf_fn.
# In closure language, the driver calls a function factory that takes T and
# returns x -> lnf_typed_fn(T, which_f, x). That returned function has x as
# its only explicit argument, but it still has access to T and which_f from
# the scope in which the closure was created.
stieltjes_kwargs = (;)

if @isdefined offset
    stieltjes_kwargs = (; stieltjes_kwargs..., offset=offset)
end

if @isdefined j_max
    stieltjes_kwargs = (; stieltjes_kwargs..., j_max=j_max)
end

if @isdefined epsilon
    stieltjes_kwargs = (; stieltjes_kwargs..., epsilon=epsilon)
end

stieltjes_final_user_fn = (n) ->
CustomGaussQuadrature.stieltjes_a_vec_b_vec_final_fn(n, μ₀_input, lnf_typed_fn, which_f, a, b; stieltjes_kwargs...)

@assert n ≥ 1
if n == 1
    a_vec_stieltjes, b_vec_stieltjes, nbits_stieltjes, r = 
    stieltjes_final_user_fn(2);
    nodes_stieltjes = convert(Double64, a_vec_stieltjes[1])
    weights_stieltjes = μ₀_double64
else
    a_vec_stieltjes, b_vec_stieltjes, nbits_stieltjes, r = 
    stieltjes_final_user_fn(n);

    nodes_stieltjes, weights_stieltjes = 
    stieltjes_custom_gauss_quad_all_fn(n, μ₀_double64, a_vec_stieltjes, b_vec_stieltjes, a, b);
end

stieltjes_nodes = nodes_stieltjes
stieltjes_weights = weights_stieltjes
stieltjes_a_vec = a_vec_stieltjes
stieltjes_b_vec = b_vec_stieltjes
stieltjes_nbits = nbits_stieltjes