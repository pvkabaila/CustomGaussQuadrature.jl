# Test the code from my CustomGaussQuadrature package
# for computing the Gauss rule using moment determinants.

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
# Test of gauss_quad_moment_dets_scr.jl
#****************************************************
###################################################
# Use
# moment_fn = moment_stored_fn
# where moment_stored_fn has been exported from
# the module CustomGaussQuadrature
###################################################

# Try which_f = ["scaled chi pdf", [0,Inf], m]
# Julia's BigFloat number type corresponds to 
# 256 bits i.e.
# nbits <- 256
# in my R code.
println("Test of the version of the CustomGaussQuadrature package")
println("on my computer, not the version in the Julia General Registry.")
println(" ")

println("✔ means exact agreement with R computed results", "\n")

moment_fn = moment_stored_fn

T = BigFloat;
precision(T)
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
# μ_offsetvec_fn(BigFloat, which_f, n) matches the result of my R code  ✔

n = 10;
println("n = ", n, "\n")
μ_offsetvec = μ_offsetvec_fn(BigFloat, moment_fn, which_f, n);
println("Δ_fn(μ_offsetvec, n, n) = ", "\n", Δ_fn(μ_offsetvec, n, n), "\n")
# Δ_fn(μ_offsetvec, n, n) matches the results of my R code  ✔

println("Δ′_fn(μ_offsetvec, n, n) = ", "\n", Δ′_fn(μ_offsetvec, n, n), "\n")
# Δ′_fn(μ_offsetvec, n, n) matches the results of my R code ✔

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
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL] 

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
println("max_abs_rel_error_nodes = ", max_abs_rel_error_nodes, "\n")

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

println("✔ means exact agreement with R computed results")

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

which_f = ["Legendre", [-1, 1]];
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

which_f = ["Generalized Laguerre", [0, Inf], 3]::Vector{Any};
println("which_f = ", which_f)
α = 3.1;
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
end

moment_fn = moment_weibull_pdf_fn


k = 2.0
which_f = ["weibull pdf", [0, Inf], k]
k = which_f[3]
println("k = ", k)

T = Float64
r = 4
moment_fn(T, which_f, r)

T = BigFloat
moment_fn(T, which_f, r)

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

β_vec = β_vec_fn(Δ_offsetvec, n)
println("β_vec = ")
for element in β_vec
    println(element)
end

#---------------------------------------------------------------
n= 33;
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

