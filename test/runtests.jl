using CustomGaussQuadrature
using Test


@testset "CustomGaussQuadrature.jl" begin

#****************************************************
# Test of gauss_quad_moment_dets_scr.jl
#****************************************************

    T = BigFloat;
    m = 160;
    which_f = ["scaled chi pdf", [0,Inf], m];

    n = 4;
    μ_offsetvec = μ_offsetvec_fn(BigFloat, which_f, n);

    Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
    Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);

    α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
    β_vec = β_vec_fn(Δ_offsetvec, n);

#----------------------------------------------------
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL] 

T = BigFloat;
which_f = ["Generalized Laguerre", [0, Inf], 1]; 
n = 4;

μ_offsetvec = μ_offsetvec_fn(T, which_f, n);

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

nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(which_f, n);
nodes = convert(Vector{Float64}, nodes_Double64);
weights = convert(Vector{Float64}, weights_Double64);


setprecision(BigFloat, 256, base=2);
n= 4;
which_f = ["Hermite", [-Inf, Inf]]
nodes, weights = custom_gauss_quad_all_fn(which_f, n);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# chemistry example for n = 4
which_f = ["chemistry example", [0, Inf]]::Vector{Any};
n = 4;
nodes, weights = custom_gauss_quad_all_fn(which_f, n);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# Generalized Laguerre
n= 4;
α = convert(BigFloat, 3);

which_f = ["Generalized Laguerre", [0, Inf], 3]::Vector{Any};
nodes, weights = custom_gauss_quad_all_fn(which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=false and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(which_f, n, false, true);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=false

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes, weights = custom_gauss_quad_all_fn(which_f, n, true);

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
n= 4;
nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(which_f, n, true, true);


#****************************************************
#  Test of gauss_quad_stjieltjes_scr.jl
#****************************************************

m = 160;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);
n = 5;

stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

a_vec, b_vec, μ₀, nbits = 
step1_fn(which_f, n);

# nodes and weights

stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

nodes, weights = 
custom_gauss_quad_all_fn(which_f, n);

#------------------------------

upto_n = true;

T = BigFloat;
a = convert(T,0);
b = Inf;

stjieltjes_nodes_upto_n, stjieltjes_weights_upto_n = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b, upto_n);

nodes_upto_n, weights_upto_n = 
custom_gauss_quad_all_fn(which_f, n, upto_n);

#--------------------------------------------------
# chemistry example with n = 5

n = 5;

which_f = ["chemistry example", [0, Inf]];
nodes, weights = custom_gauss_quad_all_fn(which_f, n);

lnf_fn = x -> lnf_chemistry_fn(T, x);
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

stjieltjes_nodes_Float64 = convert(Vector{Float64}, stjieltjes_nodes);
stjieltjes_weights_Float64 = convert(Vector{Float64}, stjieltjes_weights);

#------------------------------------------------
# Test the function stjieltjes_custom_gauss_quad_all_fn
# Hermite
setprecision(BigFloat, 256, base=2);
n = 4;
which_f = ["Hermite", [-Inf, Inf]];
lnf_fn = lnf_hermite_fn;
a = -Inf;
b = Inf;
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, k = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

#------------------------------------------------
# Test the function stjieltjes_custom_gauss_quad_all_fn
# Generalized Laguerre

setprecision(BigFloat, 256, base=2);
n= 4;

T = BigFloat;
α = convert(T, 1);
which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
a = convert(T,0);
b = Inf;
lnf_fn = x -> lnf_laguerre_fn(x, α);

μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1); 
stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, k = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

end
