using CustomGaussQuadrature
using CustomGaussQuadrature: Double64, μ_offsetvec_fn
using SpecialFunctions
using Test

@testset "String spec support" begin

@test materialize_scalar_spec_fn(Double64, "1.1") == parse(Double64, "1.1")
@test materialize_scalar_spec_fn(Double64, "1.1") != convert(Double64, 1.1)
@test isinf(materialize_scalar_spec_fn(Double64, BigFloat(Inf)))
@test materialize_scalar_spec_fn(Double64, BigFloat(Inf)) > zero(Double64)
@test isinf(materialize_scalar_spec_fn(Double64, BigFloat(-Inf)))
@test materialize_scalar_spec_fn(Double64, BigFloat(-Inf)) < zero(Double64)

numeric_gll = ["generalized laguerre", [0, Inf], 1.0]
string_gll = ["generalized laguerre", ["0", "Inf"], "1"]

numeric_nodes, numeric_weights = custom_gauss_quad_all_fn(moment_stored_fn, numeric_gll, 4)
string_nodes, string_weights = custom_gauss_quad_all_fn(moment_stored_fn, string_gll, 4)

@test all(isapprox.(string_nodes, numeric_nodes; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(string_weights, numeric_weights; rtol=1.0e-24, atol=1.0e-28))

function moment_weibull_spec_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    T_k = materialize_scalar_spec_fn(T, which_f[3])
    @assert T_k > zero(T_k)
    @assert r ≥ 0
    if r == 0
        return(one(T))
    end
    gamma(one(T) + convert(T, r) / T_k)
end

numeric_weibull = ["weibull pdf", [0, Inf], 2.0]
string_weibull = ["weibull pdf", ["0", "Inf"], "2"]

numeric_nodes, numeric_weights = custom_gauss_quad_all_fn(moment_weibull_spec_fn, numeric_weibull, 4)
string_nodes, string_weights = custom_gauss_quad_all_fn(moment_weibull_spec_fn, string_weibull, 4)

@test all(isapprox.(string_nodes, numeric_nodes; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(string_weights, numeric_weights; rtol=1.0e-24, atol=1.0e-28))

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))

@eval Main begin
    which_f = ["generalized laguerre", [0, Inf], 1.0]
    n = 4
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

stored_numeric_nodes = copy(Main.stieltjes_nodes)
stored_numeric_weights = copy(Main.stieltjes_weights)

@eval Main begin
    which_f = ["generalized laguerre", ["0", "Inf"], "1.0"]
    n = 4
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

@test all(isapprox.(Main.stieltjes_nodes, stored_numeric_nodes; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(Main.stieltjes_weights, stored_numeric_weights; rtol=1.0e-24, atol=1.0e-28))

@eval Main begin
    function lnf_user_weibull_spec_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
        @assert which_f[1] == "weibull pdf"
        T_k = materialize_scalar_spec_fn(T, which_f[3])
        @assert T_k > zero(T_k)
        @assert x > zero(T)
        log(T_k) + (T_k - one(T)) * log(x) - x^T_k
    end

    which_f = ["weibull pdf", [0, Inf], 1.0]
    n = 3
    lnf_typed_fn = lnf_user_weibull_spec_fn
    mu0 = 1
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

user_numeric_nodes = copy(Main.stieltjes_nodes)
user_numeric_weights = copy(Main.stieltjes_weights)

@eval Main begin
    which_f = ["weibull pdf", ["0", "Inf"], "1"]
    n = 3
    lnf_typed_fn = lnf_user_weibull_spec_fn
    mu0 = "1"
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

@test all(isapprox.(Main.stieltjes_nodes, user_numeric_nodes; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(Main.stieltjes_weights, user_numeric_weights; rtol=1.0e-24, atol=1.0e-28))

@eval Main begin
    which_f = ["weibull pdf", ["0", "Inf"], "1"]
    n = 3
    lnf_typed_fn = lnf_user_weibull_spec_fn
    mu0 = 1
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

@test all(isapprox.(Main.stieltjes_nodes, user_numeric_nodes; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(Main.stieltjes_weights, user_numeric_weights; rtol=1.0e-24, atol=1.0e-28))

moment_fn = moment_stored_fn
which_f = ["scaled chi pdf", [0, Inf], 160]
n = 5
lnf_typed_fn = (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3])
μ₀, μ₁ = μ_offsetvec_fn(BigFloat, moment_fn, which_f, 1)

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r =
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, which_f[2][1], which_f[2][2])

stieltjes_nodes_upto_n, stieltjes_weights_upto_n =
    stieltjes_custom_gauss_quad_all_fn(
        n,
        μ₀,
        stieltjes_a_vec,
        stieltjes_b_vec,
        BigFloat(0),
        Inf,
        true,
    )

@test length(stieltjes_nodes_upto_n) == n
@test length(stieltjes_weights_upto_n) == n
@test stieltjes_nodes_upto_n[1] isa Double64
@test isfinite(stieltjes_nodes_upto_n[end][end])
@test isfinite(stieltjes_weights_upto_n[end][end])

end