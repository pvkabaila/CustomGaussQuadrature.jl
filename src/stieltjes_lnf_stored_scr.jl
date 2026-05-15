# Stored Stieltjes driver for built-in weight functions.
#
# This file is intentionally not included from
# CustomGaussQuadrature.jl. It is meant to be run explicitly with
# include(...), typically from the REPL, tests, or an interface script.
#
# Required inputs in the caller:
#   which_f
#   n
# Optional inputs:
#   offset, j_max, epsilon
#
# Outputs in the caller:
#   stieltjes_nodes
#   stieltjes_weights
#   stieltjes_a_vec
#   stieltjes_b_vec
#   stieltjes_nbits
#   r

using CustomGaussQuadrature: Double64

a, b = which_f[2];
# Keep the support endpoints in their input representation.
# The support-quadrature helpers convert a and b to the local working type
# chosen by the underlying arithmetic path when they need endpoint values.

stored_weight_name = CustomGaussQuadrature.normalize_stored_weight_name_fn(which_f[1])

if stored_weight_name == "scaled chi pdf"
    m = materialize_integer_spec_fn(which_f[3]);
    @assert m > 0
    lnf_typed_fn = (T, which_f, x) ->
    ln_scaled_chi_pdf_fn(T, x, materialize_integer_spec_fn(which_f[3]));
    μ₀_input = T -> one(T)
elseif stored_weight_name == "hermite"
    lnf_typed_fn = (T, which_f, x) -> lnf_hermite_fn(T, x);
    μ₀_input = T -> sqrt(convert(T, π))
elseif stored_weight_name == "generalized laguerre"
    lnf_typed_fn = (T, which_f, x) -> lnf_laguerre_fn(T, x, which_f[3]);
    μ₀_input = (T, which_f) -> begin
        T_α_gl = materialize_scalar_spec_fn(T, which_f[3])
        @assert T_α_gl > -one(T_α_gl)
        gamma(T_α_gl + one(T_α_gl))
    end
elseif stored_weight_name == "generalized normal"
    @assert length(which_f[3]) == 2
    α, β = which_f[3]
    lnf_typed_fn = (T, which_f, x) -> lnf_generalized_normal_fn(T, x, α, β)
    μ₀_input = T -> one(T)
elseif stored_weight_name == "chemistry example"
    lnf_typed_fn = (T, which_f, x) -> lnf_chemistry_fn(T, x);
    μ₀_input = T -> begin
        T_1 = one(T)
        T_2 = convert(T, 2)
        T_3 = convert(T, 3)
        T_3^(-T_2 / T_3) * gamma(T_1 / T_3)
    end
else
    # DomainError means that the argument to a function
    # or constructor does not lie in the valid domain
    throw(DomainError(which_f, "invalid argument"))
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
stieltjes_a_vec_b_vec_final_fn(n, μ₀_input, lnf_typed_fn, which_f, a, b; stieltjes_kwargs...)

@assert n ≥ 1
if n == 1
    stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_final_user_fn(2);
    stieltjes_nodes = convert(Double64, stieltjes_a_vec[1])
    stieltjes_weights = μ₀_double64
else
    stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_final_user_fn(n);

    stieltjes_nodes, stieltjes_weights = 
    stieltjes_custom_gauss_quad_all_fn(n, μ₀_double64, stieltjes_a_vec, stieltjes_b_vec, a, b);
end
