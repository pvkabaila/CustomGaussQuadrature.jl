# test_CustomGaussQuad_scr.jl
# The purpose of this script is to test the
# module CustomGaussQuad.jl, which is a 
# translated version of my R code, with a few
# enhancements

using .CustomGaussQuad
using GaussQuadrature
using SpecialFunctions
using OffsetArrays
using LinearAlgebra
using GenericLinearAlgebra
using DoubleFloats
using Printf

#----------------------------------------------------
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

@time β_vec = β_vec_fn(μ_offsetvec, Δ_offsetvec, n)
# β_vec_fn(μ_offsetvec, Δ_offsetvec, n) matches the results of my R code  ✔


#----------------------------------------------------
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL] 

T = BigFloat
which_f = ["Generalized Laguerre", [0, Inf], 1]  
n = 10

μ_offsetvec = μ_offsetvec_fn(T, which_f, n)

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n)
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n)

α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n)
β_vec =  β_vec_fn(μ_offsetvec, Δ_offsetvec, n)

a, b = laguerre_coefs(n, convert(T, 1))

# The parent function converts an Offset Array into an ordinary Array
parent(α_offsetvec) - a

# The parent function converts an Offset Array into an ordinary Array
sqrt.(β_vec)- b[2:n]

#----------------------------------------------

# Default precision of BigFloat is 256 bits
precision(BigFloat)
2.0^(-precision(BigFloat))

# Store default precision of BigFloat in old_precision
old_precision = precision(BigFloat)

# Set precision of BigFloat to 50 decimal digits
setprecision(BigFloat, 50, base=10)
precision(BigFloat)
2.0^(-precision(BigFloat))

# Restore precision of BigFloat to the default
setprecision(BigFloat, old_precision, base=2) 
precision(BigFloat)

#--------------------------------------


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
@time nodes_BigFloat132, weights_BigFloat132 = custom_gauss_quad(BigFloat, which_f, n);

abs_error_nodes = convert(Vector{Float64}, abs.(nodes_Double64 - nodes_BigFloat132));
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", max_abs_error_nodes)

abs_rel_error_nodes = convert(Vector{Float64}, (abs.((nodes_Double64 - nodes_BigFloat132) ./ nodes_BigFloat132)));
max_abs_rel_error_nodes = maximum(abs_rel_error_nodes);
println("max_abs_rel_error_nodes = ", max_abs_rel_error_nodes)

abs_rel_error_weights = convert(Vector{Float64}, (abs.((weights_Double64 - weights_BigFloat132) ./ weights_BigFloat132)));
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights)

row_number_last = ceil(Int64, length(nodes_Float64)/3);
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
which_f = ["scaled chi pdf", [0,Inf], 2]::Vector{Any}
n= 33;
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
    @printf "%.18f " nodes[3*(row_number-1) + 1]
    if 3*(row_number-1) + 2 ≤ n
        @printf "%.18f " nodes[3*(row_number-1) + 2]
    end
    if 3*(row_number-1) + 3 ≤ n
        @printf "%.18f \n" nodes[3*(row_number-1) + 3]
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
n= 100;
println("n = ", n)
@time nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
which_f = ["Hermite", [-Inf, Inf]]
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n);

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))




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

#------------------------------------------------


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
n= 100;
println("n = ", n)
@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);
which_f = ["Legendre", [-1, 1]]
@time nodes, weights = custom_gauss_quad_all_fn(which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))



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

#-------------------------------------------------------










































