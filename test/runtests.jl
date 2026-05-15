using CustomGaussQuadrature
using CustomGaussQuadrature: μ_offsetvec_fn, Δ_offsetvec_fn, Δ′_offsetvec_fn, α_offsetvec_fn, β_vec_fn, step1_fn
using Test
using SpecialFunctions

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))


@testset "CustomGaussQuadrature.jl" begin

#****************************************************
# Test of gauss_quad_moment_dets_scr.jl
#****************************************************

moment_fn = moment_stored_fn

T = BigFloat;
m = 160;
which_f = ["scaled chi pdf", [0,Inf], m]::Vector{Any};

n = 4;
μ_offsetvec = μ_offsetvec_fn(BigFloat, moment_fn, which_f, n);
    
Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);

α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
β_vec = β_vec_fn(Δ_offsetvec, n);

#----------------------------------------------------
# Try which_f = ["generalized laguerre", [0, Inf], α_gl] 

T = BigFloat;
which_f = ["generalized laguerre", [0, Inf], 1]; 
n = 4;

μ_offsetvec = μ_offsetvec_fn(T, moment_fn, which_f, n);

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);

α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
β_vec =  β_vec_fn(Δ_offsetvec, n);

parent(α_offsetvec);
sqrt.(β_vec);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# scaled chi pdf weight function 

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n = 4;

nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{Float64}, nodes_Double64);
weights = convert(Vector{Float64}, weights_Double64);

# Test the function custom_gauss_quad_all_fn with 
# Hermite weight function 

setprecision(BigFloat, 256, base=2);
n= 4;
which_f = ["hermite", [-Inf, Inf]]
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# chemistry example for n = 4
which_f = ["chemistry example", [0, Inf]]::Vector{Any};
n = 4;
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# Generalized Laguerre
n= 4;
α_spec = "3.1";
α = parse(BigFloat, α_spec);

which_f = ["generalized laguerre", [0, Inf], α_spec]::Vector{Any};
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n, true);

#----------------------------------------------
# Test using moment_fn = moment_weibull_pdf_fn

function moment_weibull_pdf_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
@assert which_f[1] == "weibull pdf"
T_k = materialize_scalar_spec_fn(T, which_f[3])
@assert T_k > convert(T, 0)
@assert r ≥ 0
if r == 0
    return(convert(T, 1))
end
T_1 = convert(T, 1)
T_r = convert(T, r)
gamma(T_1 + (T_r/T_k))
end

moment_fn = moment_weibull_pdf_fn

which_f = ["weibull pdf", [0, Inf], "5.1"]

Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);

#----------------------------------------------
# Test a user-defined weight with a 2-parameter vector in which_f[3]

function moment_gamma_pdf_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
@assert which_f[1] == "gamma pdf"
α_spec, θ_spec = which_f[3]
T_α = materialize_scalar_spec_fn(T, α_spec)
T_θ = materialize_scalar_spec_fn(T, θ_spec)
@assert T_α > zero(T)
@assert T_θ > zero(T)
@assert r ≥ 0
T_r = convert(T, r)
T_θ^T_r * gamma(T_α + T_r) / gamma(T_α)
end

moment_fn = moment_gamma_pdf_fn
which_f = ["gamma pdf", [0, Inf], ["2.5", "1.7"]]
n = 4

nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n)
@test length(nodes) == n
@test length(weights) == n
@test all(nodes .>= 0)
@test all(weights .> 0)

#****************************************************
#  Test of gauss_quad_stieltjes_scr.jl
#****************************************************

m = 160;
T = BigFloat;
which_f = ["scaled chi pdf", [0,Inf], m];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3]);

moment_fn = moment_stored_fn

μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
n = 5;

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);

a_vec, b_vec, μ₀, nbits = 
step1_fn(moment_fn, which_f, n);

# nodes and weights

stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b);  

nodes, weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n);

#------------------------------

upto_n = true;

T = BigFloat;
a, b = which_f[2];

stieltjes_nodes_upto_n, stieltjes_weights_upto_n = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b, upto_n);  

nodes_upto_n, weights_upto_n = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, upto_n);

#--------------------------------------------------
# chemistry example with n = 5

n = 5;

which_f = ["chemistry example", [0, Inf]];
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_chemistry_fn(T, x);
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b);  

stieltjes_nodes_Float64 = convert(Vector{Float64}, stieltjes_nodes);
stieltjes_weights_Float64 = convert(Vector{Float64}, stieltjes_weights);

#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Hermite
setprecision(BigFloat, 256, base=2);
n = 4;
which_f = ["hermite", [-Inf, Inf]];
T = BigFloat;
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_hermite_fn(T, x);
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);

stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b);  

#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Generalized Laguerre

setprecision(BigFloat, 256, base=2);
n= 4;

T = BigFloat;
α = convert(T, 1);
which_f = ["generalized laguerre", [0, Inf], α]::Vector{Any};
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_laguerre_fn(T, x, which_f[3]);

μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);

stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b); 

#------------------------------------------------
# Test a user-defined Stieltjes weight with a 2-parameter vector in which_f[3]

setprecision(BigFloat, 256, base=2);
n = 4;
which_f = ["gamma pdf", [0, Inf], ["2.5", "1.7"]];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> begin
    α_spec, θ_spec = which_f[3]
    T_α = materialize_scalar_spec_fn(T, α_spec)
    T_θ = materialize_scalar_spec_fn(T, θ_spec)
    (T_α - one(T)) * log(x) - x / T_θ - loggamma(T_α) - T_α * log(T_θ)
end
μ₀ = (T, which_f) -> one(T)

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀(BigFloat, which_f), stieltjes_a_vec, stieltjes_b_vec, a, b);

@test length(stieltjes_nodes) == n
@test length(stieltjes_weights) == n
@test all(stieltjes_nodes .>= 0)
@test all(stieltjes_weights .> 0)

#------------------------------------------------
# Regression test for the stored generalized normal weight function

setprecision(BigFloat, 256, base=2);
n = 2;
which_f = ["generalized normal", [-Inf, Inf], ["1.9", "2.8"]];

nodes_moment, weights_moment = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);

a, b = which_f[2];
@eval Main begin
    which_f = ["generalized normal", [-Inf, Inf], ["1.9", "2.8"]]
    n = 2
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

@test length(nodes_moment) == n
@test length(weights_moment) == n
@test length(Main.stieltjes_nodes) == n
@test length(Main.stieltjes_weights) == n
@test all(isapprox.(Main.stieltjes_nodes, nodes_moment; rtol=1.0e-15, atol=1.0e-15))
@test all(isapprox.(Main.stieltjes_weights, weights_moment; rtol=1.0e-17, atol=1.0e-17))

include("test_stieltjes_driver_owned_T_scr.jl")
include("test_string_spec_support_scr.jl")

end
