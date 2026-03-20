using CustomGaussQuadrature
using CustomGaussQuadrature: μ_offsetvec_fn, Δ_offsetvec_fn, Δ′_offsetvec_fn, α_offsetvec_fn, β_vec_fn, step1_fn
using Test
using SpecialFunctions


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
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL] 

T = BigFloat;
which_f = ["Generalized Laguerre", [0, Inf], 1]; 
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
which_f = ["Hermite", [-Inf, Inf]]
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
α = 3.1;

which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=false and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, false, true);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=false

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n, true);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, true, true);

#----------------------------------------------
# Test using moment_fn = moment_weibull_pdf_fn

function moment_weibull_pdf_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
@assert which_f[1] == "weibull pdf"
k = which_f[3]
@assert k > 0
T_k = parse(T, string(k))
@assert r ≥ 0
if r == 0
    return(convert(T, 1))
end
T_1 = convert(T, 1)
T_r = convert(T, r)
gamma(T_1 + (T_r/T_k))
end

moment_fn = moment_weibull_pdf_fn

which_f = ["weibull pdf", [0, Inf], 2.1]

Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);

#****************************************************
#  Test of gauss_quad_stieltjes_scr.jl
#****************************************************

m = 160;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];

moment_fn = moment_stored_fn

μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
n = 5;

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

a_vec, b_vec, μ₀, nbits = 
step1_fn(moment_fn, which_f, n);

# nodes and weights

stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);

nodes, weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n);

#------------------------------

upto_n = true;

T = BigFloat;
a = convert(T,0);
b = Inf;

stieltjes_nodes_upto_n, stieltjes_weights_upto_n = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b, upto_n);

nodes_upto_n, weights_upto_n = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, upto_n);

#--------------------------------------------------
# chemistry example with n = 5

n = 5;

which_f = ["chemistry example", [0, Inf]];
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

lnf_fn = x -> lnf_chemistry_fn(T, x);
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);

stieltjes_nodes_Float64 = convert(Vector{Float64}, stieltjes_nodes);
stieltjes_weights_Float64 = convert(Vector{Float64}, stieltjes_weights);

#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Hermite
setprecision(BigFloat, 256, base=2);
n = 4;
which_f = ["Hermite", [-Inf, Inf]];
lnf_fn = lnf_hermite_fn;
a = -Inf;
b = Inf;
T = BigFloat;

stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, k = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);

#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Generalized Laguerre

setprecision(BigFloat, 256, base=2);
n= 4;

T = BigFloat;
α = convert(T, 1);
which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
a = convert(T,0);
b = Inf;
lnf_fn = x -> lnf_laguerre_fn(x, α);

μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, k = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);

end
