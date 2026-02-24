# Test the code from my CustomGaussQuadrature package
# for computing the Gauss rule using the Stjieltjes
# procedure.

# The path to this package is the following:
# RESEARCH - NUMERICAL METHODS/QUADRATURE - Custom Gauss/
# CustomGaussQuadrature - Julia package/CustomGaussQuadrature

# The latest version of my Julia package CustomGaussQuadrature is in the folder:
# \RESEARCH - NUMERICAL METHODS\QUADRATURE\Custom GAUSS\CustomGaussQuadrature - Julia package\
# CustomGaussQuadrature\

# To run the code in this package, the first step in VS code is 
# File > Open Folder... > open the above folder

# Use CTRL + Shift + P to get to the Command Palette in VS code and choose
# Julia: Start REPL
# to start the Julia REPL.

# Use 
# Go to package mode by typing ]
# Activate the package environment by running the 
# activate . command
# to activate the package. Exit package mode using Ctrl + C

# Then use the following command to run this script:
# julia> include("test/test_CustomGaussQuadrature_scr.jl")

# Copy the output at the REPL into a 
# text document. Do NOT execute this 
# code line-by-line in the REPL, as this introduces
# ugly extra spaces.

# Have a look at the README.md file by opening this file in VS code and then using
# CTRL + K      V 
# to get a preview of this file.

using CustomGaussQuadrature
# This script tests the version on my computer, not the
# version in the Julia General Registry. See path.
pathof(CustomGaussQuadrature)

using Printf
using Plots
using GaussQuadrature
using SpecialFunctions
using QuadGK
using DoubleFloats

#****************************************************
# Test of gauss_quad_moment_dets_scr.jl
#****************************************************
# Try which_f = ["scaled chi pdf", [0,Inf], m]
# Julia's BigFloat number type corresponds to 
# 256 bits i.e.
# nbits <- 256
# in my R code.

println("Test of the version of the CustomGaussQuadrature package")
println("on my computer, not the version in the Julia General Registry.")
println(" ")

println("moment_fn = moment_stored_fn;")
moment_fn = moment_stored_fn;

T = BigFloat;
m = 160;
which_f = ["scaled chi pdf", [0,Inf], m]::Vector{Any};
println("T = ", T, ",    which_f = ", which_f, "\n")

n = 3;
println("n = ", n, "\n")
μ_offsetvec = μ_offsetvec_fn(BigFloat, moment_fn, which_f, n);
println("μ_offsetvec = ")
for element in μ_offsetvec
    println(element)
end
print("\n")

n = 10;
println("n = ", n, "\n")
μ_offsetvec = μ_offsetvec_fn(BigFloat, moment_fn, which_f, n);
println("Δ_fn(μ_offsetvec, n, n) = ", "\n", Δ_fn(μ_offsetvec, n, n), "\n")

println("Δ′_fn(μ_offsetvec, n, n) = ", "\n", Δ′_fn(μ_offsetvec, n, n), "\n")

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
println("Δ_offsetvec = ")
for element in Δ_offsetvec
    println(element)
end
print("\n")

Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);
println("Δ′_offsetvec = ")
for element in Δ′_offsetvec
    println(element)
end
print("\n")


α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
println("α_offsetvec = ")
for element in α_offsetvec
    println(element)
end
print("\n")
# α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n) matches the results of my R code  ✔

β_vec = β_vec_fn(Δ_offsetvec, n)
println("β_vec = ")
for element in β_vec
    println(element)
end
# β_vec_fn(Δ_offsetvec, n) matches the results of my R code  ✔

println("\n", "------------------------------------------------------")
#----------------------------------------------------
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL], where α_GGL=1

T = BigFloat;
which_f = ["Generalized Laguerre", [0, Inf], 1]::Vector{Any}; 
println("T = ", T, ",   which_f = ", which_f, "\n")

n = 10;
println("n = ", n, "\n")

μ_offsetvec = μ_offsetvec_fn(T, moment_fn, which_f, n);

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);

α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
β_vec =  β_vec_fn(Δ_offsetvec, n);

println("Use the command laguerre_coefs from the GaussQuadrature package:")
println("a, b = laguerre_coefs(n, convert(T, 1))", "\n")
a, b = laguerre_coefs(n, convert(T, 1));

println("The parent function converts an Offset Array into an ordinary Array")
println("parent(α_offsetvec) - a = ")
# The parent function converts an Offset Array into an ordinary Array
for element in (parent(α_offsetvec) - a)
    println(element)
end
print("\n")

println("sqrt.(β_vec)- b[2:n] = ")
for element in (sqrt.(β_vec)- b[2:n])
    println(element)
end

println("\n", "------------------------------------------------------")
#--------------------------------------------------------------
# Test of the function custom_gauss_quad_fn and print the nodes
# and weights in the same format as that used for the output
# from my R package custom.gauss.quad

println("Test of the function custom_gauss_quad_fn and print the nodes")
println("and weights in the same format as that used for the output")
println("from my R package custom.gauss.quad.", "\n")

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 33;
println("n = ", n, "\n")


println("@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

setprecision(BigFloat, 132, base=2);
println("BigFloat precision set to 132 bit significand")
println("@time nodes_BigFloat132, weights_BigFloat132 = custom_gauss_quad_fn(BigFloat, moment_fn, which_f, n);")
@time nodes_BigFloat132, weights_BigFloat132 = custom_gauss_quad_fn(BigFloat, moment_fn, which_f, n);
print("\n")

abs_error_nodes = convert(Vector{Float64}, abs.(nodes_Double64 - nodes_BigFloat132));
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", max_abs_error_nodes)

abs_rel_error_nodes = convert(Vector{Float64}, (abs.((nodes_Double64 - nodes_BigFloat132) ./ nodes_BigFloat132)));
max_abs_rel_error_nodes = maximum(abs_rel_error_nodes);
println("max_abs_rel_error_nodes = ", max_abs_rel_error_nodes)

abs_rel_error_weights = convert(Vector{Float64}, (abs.((weights_Double64 - weights_BigFloat132) ./ weights_BigFloat132)));
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights, "\n")

row_number_last = ceil(Int64, length(nodes_Double64)/3);
println("Double64 nodes:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 1]
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 2]
    @printf "%.17f \n" nodes_Double64[3*(row_number-1) + 3]
end
print("\n")
println("Double64 weights:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights_Double64[3*(row_number-1) + 1]
    @printf "%.16e " weights_Double64[3*(row_number-1) + 2]
    @printf "%.16e \n" weights_Double64[3*(row_number-1) + 3]
end

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# scaled chi pdf
which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n = 33;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{Float64}, nodes_Double64);
weights = convert(Vector{Float64}, weights_Double64);
# Print the nodes and weights in the same format as that used for the output
# from my R package custom.gauss.quad
print("\n")
row_number_last = ceil(Int64, length(nodes)/3);

println("Double64 nodes after conversion to Float64:")
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
print("\n")

println("Double64 weights after conversion to Float64:")
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

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Hermite
which_f = ["Hermite", [-Inf, Inf]];
println("which_f = ", which_f)

println("precision(BigFloat) = ",precision(BigFloat))
println("setprecision(BigFloat, 256, base=2) ", "\n")
setprecision(BigFloat, 256, base=2);

n= 5;
println("n = ", n, "\n")

println("Use the command hermite from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);")
@time nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

println("        nodes_BigFloat (GaussQuadrature)               weights_BigFloat (GaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

println("         nodes (CustomGaussQuadrature)                weights (CustomGaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end
print("\n")

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# chemistry example
which_f = ["chemistry example", [0, Inf]]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 15;
println("n = ", n, "\n")
println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

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
println("Print the nodes and weights in the same format as that used in")
println("Table 2.2 of Gautschi (1983):")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " nodes[i]
    @printf "%.15e  \n" weights[i]
end

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Legendre
nwhich_f = ["Legendre", [-1, 1]];
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("Use the command legendre from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);")
@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
println("precision(BigFloat) = ", precision(BigFloat))
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
print("\n")

println("        nodes_BigFloat (GaussQuadrature)               weights_BigFloat (GaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

println("         nodes (CustomGaussQuadrature)                weights (CustomGaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end
print("\n")

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Generalized Laguerre
α = 3.1;
which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
println("which_f = ", which_f)

println("α = ", α, "\n")

n= 5;
println("n = ", n, "\n")

println("Use the command laguerre from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);")
@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
print("\n")

println("        nodes_BigFloat (GaussQuadrature)               weights_BigFloat (GaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

println("         nodes (CustomGaussQuadrature)                weights (CustomGaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end
print("\n")

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=false and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights =")
println("custom_gauss_quad_all_fn(which_f, n, false, true);")
@time nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, false, true);
print("\n")

println("nodes_Double64 = ")
for element in nodes_Double64
    println(element)
end
print("\n")

println("weights_Double64 = ")
for element in weights_Double64
    println(element)
end
print("\n")

println("max_abs_error_nodes = ", max_abs_error_nodes)
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights)

println("\n", "------------------------------------------------------")
#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=false

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n, true);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n, true);
print("\n")


@printf "%2d     " 1
    @printf "node:   %.17f     " convert(BigFloat, nodes[1])
    @printf "weight: %.17f  \n" convert(BigFloat, weights[1])
    @printf "\n"
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

println("\n", "------------------------------------------------------")
#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=true

println("Test the function custom_gauss_quad_all_fn with") 
println("upto_n=true and extra_check=true", "\n")

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("@time nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights =")
println("custom_gauss_quad_all_fn(which_f, n, true, true);")
@time nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, true, true);
print("\n")


@printf "%2d     " 1
    @printf "node:   %.17f     " convert(BigFloat, nodes[1])
    @printf "weight: %.17f  \n" convert(BigFloat, weights[1])
    @printf "\n"
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

println("max_abs_error_nodes = ",  max_abs_error_nodes)
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights)


###################################################
# Use, for example,
# moment_fn = moment_weibull_pdf_fn
###################################################
println("\n", "------------------------------------------------------")
println("Test the idea using the Weibull pdf, with scale parameter 1")
println("and shape parameter k > 0, weight function")

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
end;

println("moment_fn = moment_weibull_pdf_fn;")
moment_fn = moment_weibull_pdf_fn;

k = 2.0;
which_f = ["weibull pdf", [0, Inf], k];
println("which_f = ", which_f, ",  k=", k, "\n")

#---------------------------------------------------------------
println("\n", "------------------------------------------------------")
n= 63;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64 = ")
println("custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")
row_number_last = ceil(Int64, length(nodes_Double64)/3);
println("Double64 nodes:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 1]
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 2]
    @printf "%.17f \n" nodes_Double64[3*(row_number-1) + 3]
end
print("\n")
println("Double64 weights:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights_Double64[3*(row_number-1) + 1]
    @printf "%.16e " weights_Double64[3*(row_number-1) + 2]
    @printf "%.16e \n" weights_Double64[3*(row_number-1) + 3]
end

#-------------------------------------------------
# Plot the cdf corresponding to the Gauss quadrature 
# rule and the actual cdf of a Weibull distribution 
# with scale parameter λ = 1 and shape parameter k 
# (k > 0)

# To include a plot pane within the VS code editor window,
# open the command palette with Ctrl + Shift + P
# and enter the command
# Julia: Enable Plot Pane

x_vec = nodes_Double64;
prob_vec = weights_Double64;
x_lo = 0.0;
x_hi = 2.5;
plot_cdf_discrete_rv_fn(x_vec, prob_vec, x_lo, x_hi) 
x_grid = range(x_lo, x_hi, length=200);
y_grid = weibull_cdf_fn.(x_grid, k);
plot!(x_grid, y_grid, 
title="Weibull pdf weight function with scale parameter λ=1.0 and k=$k",
 titlefont=font(10))


#****************************************************
#  Test of gauss_quad_stjieltjes_scr.jl
#****************************************************
# Final values of a_vec and b_vec

println("moment_fn = moment_stored_fn;")
moment_fn = moment_stored_fn;

m = 160;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
println("scaled chi pdf weight function, with m = ", m)
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
n = 5;
println("number of Gauss quadrature nodes n = ", n)

@time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r)

@time "step1_fn" a_vec, b_vec, μ₀, nbits = 
step1_fn(moment_fn, which_f, n);

max_abs_diff_a_vec = maximum(abs.(stjieltjes_a_vec - a_vec));
println("maximum(abs.(stjieltjes_a_vec - a_vec)) = ", convert(Float64, max_abs_diff_a_vec))

max_abs_rel_diff_b_vec = maximum(abs.((stjieltjes_b_vec - b_vec) ./ b_vec));
println("maximum(abs.((stjieltjes_b_vec - b_vec) ./ b_vec)) = ", convert(Float64, max_abs_rel_diff_b_vec))

# nodes and weights

@time "stjieltjes_custom_gauss_quad_all_fn" stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);

@time "custom_gauss_quad_all_fn" nodes, weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n);


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
custom_gauss_quad_all_fn(moment_fn, which_f, n, upto_n);

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
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

lnf_fn = x -> lnf_chemistry_fn(T, x);
println("chemistry example weight function")
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

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
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1); 

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

μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1); 
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