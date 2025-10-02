using CustomGaussQuadrature
pathof(CustomGaussQuadrature)

using Printf
using Plots
using GaussQuadrature


#****************************************************
# Test of gauss_quad_moment_dets_scr.jl
#****************************************************
# Try which_f = ["scaled chi pdf", [0,Inf], m]
# Julia's BigFloat number type corresponds to 
# 256 bits i.e.
# nbits <- 256
# in my R code.

T = BigFloat
m = 160
which_f = ["scaled chi pdf", [0,Inf], m]

n = 3
@time μ_offsetvec = μ_offsetvec_fn(BigFloat, which_f, n)
# μ_offsetvec_fn(BigFloat, which_f, n) matches the result of my R code  ✔

n = 10
@time μ_offsetvec = μ_offsetvec_fn(BigFloat, which_f, n)
@time Δ_fn(μ_offsetvec, n, n) 
# Δ_fn(μ_offsetvec, n, n) matches the results of my R code  ✔

@time Δ′_fn(μ_offsetvec, n, n) 
# Δ′_fn(μ_offsetvec, n, n) matches the results of my R code ✔

@time Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n)
@time Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n)

@time α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n)
# α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n) matches the results of my R code  ✔

@time β_vec = β_vec_fn(Δ_offsetvec, n)
# β_vec_fn(Δ_offsetvec, n) matches the results of my R code  ✔


#----------------------------------------------------
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL] 

T = BigFloat
which_f = ["Generalized Laguerre", [0, Inf], 1]  
n = 10

μ_offsetvec = μ_offsetvec_fn(T, which_f, n)

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n)
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n)

α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n)
β_vec =  β_vec_fn(Δ_offsetvec, n)

a, b = laguerre_coefs(n, convert(T, 1))

# The parent function converts an Offset Array into an ordinary Array
parent(α_offsetvec) - a

# The parent function converts an Offset Array into an ordinary Array
sqrt.(β_vec)- b[2:n]


#--------------------------------------------------------------
# Test of the function custom_gauss_quad_fn and print the nodes
# and weights in the same format as that used for the output
# from my R package custom.gauss.quad

which_f = ["scaled chi pdf", [0,Inf], 160]
n= 33;
println("n = ", n)

@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(which_f, n);

setprecision(BigFloat, 132, base=2);
println("T = BigFloat with 132 bit significand")
@time nodes_BigFloat132, weights_BigFloat132 = custom_gauss_quad_fn(BigFloat, which_f, n);

abs_error_nodes = convert(Vector{Float64}, abs.(nodes_Double64 - nodes_BigFloat132));
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", max_abs_error_nodes)

abs_rel_error_nodes = convert(Vector{Float64}, (abs.((nodes_Double64 - nodes_BigFloat132) ./ nodes_BigFloat132)));
max_abs_rel_error_nodes = maximum(abs_rel_error_nodes);
println("max_abs_rel_error_nodes = ", max_abs_rel_error_nodes)

abs_rel_error_weights = convert(Vector{Float64}, (abs.((weights_Double64 - weights_BigFloat132) ./ weights_BigFloat132)));
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights)

row_number_last = ceil(Int64, length(nodes_Double64)/3);
@printf("Double64 nodes:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 1]
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 2]
    @printf "%.17f \n" nodes_Double64[3*(row_number-1) + 3]
end
@printf("Double64 weights:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights_Double64[3*(row_number-1) + 1]
    @printf "%.16e " weights_Double64[3*(row_number-1) + 2]
    @printf "%.16e \n" weights_Double64[3*(row_number-1) + 3]
end


#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# scaled chi pdf
which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any}
n = 33;
println("n = ", n)
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(which_f, n);
nodes = convert(Vector{Float64}, nodes_Double64);
weights = convert(Vector{Float64}, weights_Double64);
# Print the nodes and weights in the same format as that used for the output
# from my R package custom.gauss.quad

row_number_last = ceil(Int64, length(nodes)/3);

println("✔ means exact agreement with R computed results")

@printf("Double64 nodes after conversion to Float64:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes[3*(row_number-1) + 1]
    if 3*(row_number-1) + 2 ≤ n
        @printf "%.17f " nodes[3*(row_number-1) + 2]
    end
    if 3*(row_number-1) + 3 ≤ n
        @printf "%.17f \n" nodes[3*(row_number-1) + 3]
    end
end
if 3*(row_number_last-1) + 3 !== n
    println("")
end

@printf("Double64 weights after conversion to Float64:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights[3*(row_number-1) + 1]
    if 3*(row_number-1) + 2 ≤ n
        @printf "%.16e " weights[3*(row_number-1) + 2]
    end
    if 3*(row_number-1) + 3 ≤ n
        @printf "%.16e \n" weights[3*(row_number-1) + 3]
    end
end
if 3*(row_number_last-1) + 3 !== n
    println("")
end


#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Hermite
precision(BigFloat)
setprecision(BigFloat, 256, base=2) 
n= 5;
println("n = ", n)
@time nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
which_f = ["Hermite", [-Inf, Inf]]
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n);
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
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
@printf "                   nodes                                  weights"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end


error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))




#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# chemistry example
which_f = ["chemistry example", [0, Inf]]::Vector{Any}
n= 15;
println("n = ", n)
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n);

# println("✔ means exact agreement with R computed results:")
println("✔ means exact agreement with Table 2.2 of Gautschi (1983):")

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

#---------------------------------------------------------------------------
# Print the nodes and weights in the same format as that used in Table 2.2
# of Gautschi (1983)
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " nodes[i]
    @printf "%.15e  \n" weights[i]
end

#-----------------------------------------------------------------------------

#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Legendre
n= 5;
println("n = ", n)
@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);
which_f = ["Legendre", [-1, 1]]
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n);
precision(BigFloat)
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

@printf "                nodes_BigFloat                         weights_BigFloat"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end


@printf "                   nodes                                  weights"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

#------------------------------------------------


#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Generalized Laguerre
n= 5;
println("n = ", n)
α = convert(BigFloat, 3);
println("α = ", α)
@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);

which_f = ["Generalized Laguerre", [0, Inf], 3]::Vector{Any}  
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

@printf "                nodes_BigFloat                         weights_BigFloat"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end


@printf "                   nodes                                  weights"
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=false and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any}
n= 5;
println("n = ", n)
@time nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(which_f, n, false, true);
nodes_Double64
weights_Double64
max_abs_error_nodes
max_abs_rel_error_weights


#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=false

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any}
n= 5;
println("n = ", n)
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n, true);


@printf "%2d     " 1
    @printf "node:   %.17f     " convert(BigFloat, nodes[1])
    @printf "weight: %.17f  \n" convert(BigFloat, weights[1])
for i in 2:n
    @printf "%2d     \n" i
    @printf "nodes:   "
    for j in 1:i
        @printf "%.17f     " convert(BigFloat, nodes[i][j])
    end
    @printf "\n"
    @printf "weights: "
    for j in 1:i
        @printf "%.17f    " convert(BigFloat, weights[i][j])
    end
    @printf "\n \n"
end


#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any}
n= 5;
println("n = ", n)
@time nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(which_f, n, true, true);


@printf "%2d     " 1
    @printf "node:   %.17f     " convert(BigFloat, nodes[1])
    @printf "weight: %.17f  \n" convert(BigFloat, weights[1])
for i in 2:n
    @printf "%2d     \n" i
    @printf "nodes:   "
    for j in 1:i
        @printf "%.17f     " convert(BigFloat, nodes[i][j])
    end
    @printf "\n"
    @printf "weights: "
    for j in 1:i
        @printf "%.17f    " convert(BigFloat, weights[i][j])
    end
    @printf "\n \n"
end

max_abs_error_nodes 

max_abs_rel_error_weights


#****************************************************
#  Test of gauss_quad_stjieltjes_scr.jl
#****************************************************
# Final values of a_vec and b_vec

m = 160;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
println("scaled chi pdf weight function, with m = ", m)
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);
n = 5;
println("number of Gauss quadrature nodes n = ", n)

@time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r)

@time "step1_fn" a_vec, b_vec, μ₀, nbits = 
step1_fn(which_f, n);


max_abs_diff_a_vec = maximum(abs.(stjieltjes_a_vec - a_vec));
println("maximum(abs.(stjieltjes_a_vec - a_vec)) = ", convert(Float64, max_abs_diff_a_vec))

max_abs_rel_diff_b_vec = maximum(abs.((stjieltjes_b_vec - b_vec) ./ b_vec));
println("maximum(abs.((stjieltjes_b_vec - b_vec) ./ b_vec)) = ", convert(Float64, max_abs_rel_diff_b_vec))

# nodes and weights

@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

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
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b, upto_n);

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

which_f = ["chemistry example", [0, Inf]];
nodes, weights = custom_gauss_quad_all_fn(which_f, n);

lnf_fn = x -> lnf_chemistry_fn(T, x);
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