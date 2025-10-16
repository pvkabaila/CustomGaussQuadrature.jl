# Test the code from my CustomGaussQuadrature package
# for computing the Gauss rule using the Stjieltjes
# procedure.

# The path to this package is the following:
# RESEARCH - NUMERICAL METHODS/QUADRATURE - Custom Gauss/
# CustomGaussQuadrature - Julia package/CustomGaussQuadrature

# Make this file active and execute it in REPL.
# Then copy the output at the REPL into a 
# Notepad++ text document. Do NOT execute this 
# code line-by-line in the REPL, as this introduces
# ugly extra spaces.

using CustomGaussQuadrature
# This script tests the version on my computer, not the
# version in the Julia General Registry. See path.
pathof(CustomGaussQuadrature)

using Printf
using Plots
using GaussQuadrature
using SpecialFunctions


#****************************************************
#  Test of gauss_quad_stjieltjes_scr.jl
#****************************************************
# Final values of a_vec and b_vec

m = 160;
println("scaled chi pdf weight function, with m = ", m)
which_f = ["scaled chi pdf", [0,Inf], m];
moment_fn = moment_stored_fn
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
n = 5;
println("number of Gauss quadrature nodes n = ", n)

lnf_fn = lnf_stored_fn

@time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, which_f);
println("r = ", r)

@time "step1_fn" a_vec, b_vec, μ₀, nbits = 
step1_fn(which_f, n);


max_abs_diff_a_vec = maximum(abs.(stjieltjes_a_vec - a_vec));
println("maximum(abs.(stjieltjes_a_vec - a_vec)) = ", convert(Float64, max_abs_diff_a_vec))

max_abs_rel_diff_b_vec = maximum(abs.((stjieltjes_b_vec - b_vec) ./ b_vec));
println("maximum(abs.((stjieltjes_b_vec - b_vec) ./ b_vec)) = ", convert(Float64, max_abs_rel_diff_b_vec))

# nodes and weights

@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, which_f);

@time "custom_gauss_quad_all_fn" nodes, weights = 
custom_gauss_quad_all_fn(which_f, n);

println("maximum(abs.(stjieltjes_nodes - nodes) = ", 
convert(Float64, maximum(abs.(stjieltjes_nodes - nodes))))
println("maximum(abs.((stjieltjes_weights - weights) ./ weights) = ", 
convert(Float64, maximum(abs.((stjieltjes_weights - weights) ./ weights))))

stjieltjes_nodes = convert(Vector{Float64}, stjieltjes_nodes);
stjieltjes_weights = convert(Vector{Float64}, stjieltjes_weights);
# println(" ")
@printf "           stjieltjes_nodes             stjieltjes_weights"
for i in 1:lastindex(stjieltjes_nodes)
    @printf "%2d     " i
    @printf "%.16e     " stjieltjes_nodes[i]
    @printf "%.16e  \n" stjieltjes_weights[i]
end

#------------------------------

upto_n = true;
println("upto_n = ", upto_n)

T = BigFloat;
a = convert(T,0);
b = Inf;

@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes_upto_n, stjieltjes_weights_upto_n = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, which_f, upto_n);

@time "custom_gauss_quad_all_fn" nodes_upto_n, weights_upto_n = 
custom_gauss_quad_all_fn(which_f, n, upto_n);

max_abs_diff_nodes = zeros(n);
for q in 1:n
    max_abs_diff_nodes[q] = maximum(abs.(stjieltjes_nodes_upto_n[q] - nodes_upto_n[q]));
end
println("maximum(max_abs_diff_nodes) = ", maximum(max_abs_diff_nodes))

max_abs_reldiff_weights = zeros(n);
for q in 1:n
    max_abs_reldiff_weights[q] = 
    maximum(abs.((stjieltjes_nodes_upto_n[q] - nodes_upto_n[q]) ./ nodes_upto_n[q]));
end
println("maximum(max_abs_reldiff_weights) = ", maximum(max_abs_reldiff_weights))

#---------------------------------------------------------
# Conversion of nodes_upto_n into a vector of Float64 vectors

nodes_upto_n_Float64 = Vector{Vector}(undef, n);
nodes_upto_n_Float64[1] = [convert(Float64, nodes_upto_n[1])];
for q in 2:n
    nodes_upto_n_Float64[q] = convert(Vector{Float64}, nodes_upto_n[q]);
end
println("nodes_upto_n_Float64")
nodes_upto_n_Float64

#---------------------------------------------------------
# Conversion of stjieltjes_nodes_upto_n into a vector of Float64 vectors

stjieltjes_nodes_upto_n_Float64 = Vector{Vector}(undef, n);
stjieltjes_nodes_upto_n_Float64[1] = [convert(Float64, stjieltjes_nodes_upto_n[1])];
for q in 2:n
    stjieltjes_nodes_upto_n_Float64[q] = convert(Vector{Float64}, stjieltjes_nodes_upto_n[q]);
end
println("stjieltjes_nodes_upto_n_Float64")
stjieltjes_nodes_upto_n_Float64

@printf "stjieltjes_nodes_upto_n_Float64 .- nodes_upto_n_Float64"
show(stdout, "text/plain", stjieltjes_nodes_upto_n_Float64 .- nodes_upto_n_Float64)

#---------------------------------------------------------
# Conversion of weights_upto_n into a vector of Float64 vectors

weights_upto_n_Float64 = Vector{Vector}(undef, n);
weights_upto_n_Float64[1] = [convert(Float64, weights_upto_n[1])];
for q in 2:n
    weights_upto_n_Float64[q] = convert(Vector{Float64}, weights_upto_n[q]);
end
println("weights_upto_n_Float64")
weights_upto_n_Float64


#---------------------------------------------------------
# Conversion of stjieltjes_weights_upto_n into a vector of Float64 vectors

stjieltjes_weights_upto_n_Float64 = Vector{Vector}(undef, n);
stjieltjes_weights_upto_n_Float64[1] = [convert(Float64, stjieltjes_weights_upto_n[1])];
for q in 2:n
    stjieltjes_weights_upto_n_Float64[q] = convert(Vector{Float64}, stjieltjes_weights_upto_n[q]);
end
println("stjieltjes_weights_upto_n_Float64")
stjieltjes_weights_upto_n_Float64

@printf "stjieltjes_weights_upto_n_Float64 .- weights_upto_n_Float64"
show(stdout, "text/plain", stjieltjes_weights_upto_n_Float64 .- weights_upto_n_Float64)


#--------------------------------------------------
# Convert nodes and weights to Vector{Float64} and Vector{Float64},
# respectively.
# Convert stjieltjes_nodes and stjieltjes_weights to Vector{Float64}
# and Vector{Float64}, respectively

nodes_Float64 = convert(Vector{Float64}, nodes);
weights_Float64 = convert(Vector{Float64}, weights);

stjieltjes_nodes_Float64 = convert(Vector{Float64}, nodes);
stjieltjes_weights_Float64 = convert(Vector{Float64}, weights);

show(stdout, "text/plain", stjieltjes_weights_Float64 - weights_Float64)


#--------------------------------------------------
# The weight function considered by Gautschi (1983)

n = 15;
println("number of Gauss quadrature nodes n = ", n)

lnf_fn = lnf_stored_fn
moment_fn = moment_stored_fn
which_f = ["chemistry example", [0, Inf]];
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

lnf_fn = lnf_stored_fn;

println("chemistry example weight function")
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

@time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
stjieltjes_nodes_bf = convert(Vector{BigFloat}, stjieltjes_nodes);
stjieltjes_weights_bf = convert(Vector{BigFloat}, stjieltjes_weights);

#---------------------------------------------------------------------------
# Print the nodes and weights in the same format as that used in Table 2.2
# of Gautschi (1983)
@printf "       stjieltjes_nodes_bf        stjieltjes_weights_bf"
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " stjieltjes_nodes_bf[i]
    @printf "%.15e  \n" stjieltjes_weights_bf[i]
end

#---------
stjieltjes_nodes_Float64 = convert(Vector{Float64}, stjieltjes_nodes);
stjieltjes_weights_Float64 = convert(Vector{Float64}, stjieltjes_weights);
#---------------------------------------------------------------------------
# Print the nodes and weights in the same format as that used in Table 2.2
# of Gautschi (1983)
@printf "       stjieltjes_nodes_Float64   stjieltjes_weights_Float64"
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " stjieltjes_nodes_Float64[i]
    @printf "%.15e  \n" stjieltjes_weights_Float64[i]
end

#-----------------------------------------------




#------------------------------------------------
# Test the function stjieltjes_custom_gauss_quad_all_fn
# Hermite
setprecision(BigFloat, 256, base=2);
 
n = 15;
println("number of Gauss quadrature nodes n = ", n)

which_f = ["Hermite", [-Inf, Inf]];
lnf_fn = lnf_hermite_fn;
println("Hermite weight function")
a = -Inf;
b = Inf;
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

@time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r)
@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);


nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
@printf "                nodes_BigFloat                         weights_BigFloat"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
stjieltjes_nodes_bf = convert(Vector{BigFloat}, stjieltjes_nodes);
stjieltjes_weights_bf = convert(Vector{BigFloat}, stjieltjes_weights);
@printf "             stjieltjes_nodes_bf                             stjieltjes_weights_bf"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " stjieltjes_nodes_bf[i]
    @printf "%.33e  \n" stjieltjes_weights_bf[i]
end


max_abs_error_nodes = maximum(abs.(stjieltjes_nodes_bf - nodes_BigFloat));
println("maximum(abs.(stjieltjes_nodes_bf - nodes_BigFloat)) = ", convert(Float64, max_abs_error_nodes))

max_abs_rel_error_weights = maximum(abs.((stjieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat));
println("maximum(abs.((stjieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat)) = ", 
convert(Float64, max_abs_rel_error_weights))


#------------------------------------------------
# Test the function stjieltjes_custom_gauss_quad_all_fn
# Generalized Laguerre

setprecision(BigFloat, 256, base=2);
n= 15;
println("number of Gauss quadrature nodes n = ", n)

T = BigFloat;
α = convert(T, 1);
which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
a = convert(T,0);
b = Inf;
println("Generalized Laguerre, with α = ", α)
lnf_fn = x -> lnf_laguerre_fn(x, α);

nodes_BigFloat, weights_BigFloat = laguerre(n, α);

μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1); 
@time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r)

@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

@printf "                nodes_BigFloat                         weights_BigFloat"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
stjieltjes_nodes_bf = convert(Vector{BigFloat}, stjieltjes_nodes);
stjieltjes_weights_bf = convert(Vector{BigFloat}, stjieltjes_weights);
@printf "             stjieltjes_nodes_bf                             stjieltjes_weights_bf"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " stjieltjes_nodes_bf[i]
    @printf "%.33e  \n" stjieltjes_weights_bf[i]
end

max_abs_error_nodes = maximum(abs.(stjieltjes_nodes_bf - nodes_BigFloat));
println("maximum(abs.(stjieltjes_nodes_bf - nodes_BigFloat)) = ", convert(Float64, max_abs_error_nodes))

max_abs_rel_error_weights = maximum(abs.((stjieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat));
println("maximum(abs.((stjieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat)) = ", 
convert(Float64, max_abs_rel_error_weights))

#---------------------------------------