using CustomGaussQuadrature
using Test

@testset "CustomGaussQuadrature.jl" begin
    using .CustomGaussQuadrature
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

end
