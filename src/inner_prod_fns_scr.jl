
# This script provides Julia functions to compute an approximation
# to the inner product of two functions u: ℝ → ℝ and v: ℝ → ℝ. 
# Let f: ℝ → [0,∞) denote a specified nonnegative weight function. 
# The inner product of the two functions u and v is defined to be
#
#    ∞
#    ∫ u(x) v(x) f(x) dx.                                  
#   -∞
#
# To compute this inner product let the function g = uv, so that
# this integral becomes 
# 
#    ∞
#    ∫ g(x) f(x) dx.                                       
#   -∞
#
# We also suppose that the support of the weight function f
# is an interval [a,b], -∞ ≤ a < b ≤ ∞, so that this integral
# becomes
# 
#    b
#    ∫ g(x) f(x) dx.                                      (1) 
#    a
#
# To compute this integral, we use the transformation described
# on p.94 of Gautschi (2004). In other words, we transform the  
# support interval [a,b] to the interval [-1,1] using the 
# transformation
#
#   b                1
#   ∫ g(x) f(x) dx = ∫ g(ϕ(y)) f(ϕ(y)) ϕ′(y) dy              (2)
#   a               -1
#
# where  
#
#        ϕ(y) = (1/2) (b - a) y + (1/2) (b + a) if -∞ < a < b < ∞
#        
#        ϕ(y) = b - (1 - y) / (1 + y)           if -∞ = a < b < ∞
#   
#        ϕ(y) = a + (1 + y) / (1 - y)           if -∞ < a < b = ∞
#
#        ϕ(y) = y / (1 - y²)                    if -∞ = a < b = ∞
#
# The integral on the right-hand side of (1) is approximated using Gauss
# Legendre quadrature with r nodes. 
# 
# Reference
# Gautschi, W. (2004) Orthogonal Polynomials, Computation and Approximation. 
# Oxford University Press, Oxford.



"""
ϕ = ϕ_fn(T, y, a, b) 

Returns ϕ(y) computed using type T floating point arithmetic, where  
    ϕ(y) = (1/2) (b - a) y + (1/2) (b + a) if -∞ < a < b < ∞    
    ϕ(y) = b - (1 - y) / (1 + y)           if -∞ = a < b < ∞
    ϕ(y) = a + (1 + y) / (1 - y)           if -∞ < a < b = ∞
    ϕ(y) = y / (1 - y²)                    if -∞ = a < b = ∞ 
"""
function ϕ_fn(T, y::AbstractFloat, a::AbstractFloat, b::AbstractFloat) 
    @assert a < b
    T_1 = convert(T, 1)
    @assert -T_1 ≤ y ≤ T_1
    if -Inf < a && b < Inf
        T_2 = convert(T, 2)
        out = ((b - a)/T_2) * y + (b + a) / T_2
        return(out)
    elseif -Inf == a && b < Inf
        out = b - (T_1 - y) / (T_1 + y)
        return(out)
    elseif -Inf < a && b == Inf
        out = a + (T_1 + y) / (T_1 - y)
        return(out)
    else
        out = y / (T_1 - y^2)
        return(out)
    end
end


"""
ϕ′ = ϕ′_fn(T, y::AbstractFloat, a::AbstractFloat, b::AbstractFloat)

Returns ϕ′(y) computed using type T floating point arithmetic, where
ϕ′(y) = (1/2) (b - a)               if -∞ < a < b < ∞    
ϕ′(y) = 2 / (1 + y)²                if -∞ = a < b < ∞
ϕ′(y) = 2 / (1 - y)²                if -∞ < a < b = ∞
ϕ′(y) = (1 + y²) / (1 - y²)²        if -∞ = a < b = ∞ 
"""
function ϕ′_fn(T, y::AbstractFloat, a::AbstractFloat, b::AbstractFloat) 
    @assert a < b
    T_1 = convert(T, 1)
    T_2 = convert(T, 2)
    @assert -T_1 ≤ y ≤ T_1
    if -Inf < a && b < Inf
        out = (b - a) / T_2
        return(out)
    elseif -Inf == a && b < Inf
        out =  T_2 / (T_1 + y)^2
        return(out)
    elseif -Inf < a && b == Inf
        out =  T_2 / (T_1 - y)^2
        return(out)
    else
        out = (T_1 + y^2) / (T_1 - y^2)^2
        return(out)
    end
end


"""
ln_ϕ′ = ln_ϕ′_fn(T, y::AbstractFloat, a::AbstractFloat, b::AbstractFloat)

Returns log(ϕ′(y)) computed using type T floating point arithmetic,
when it is NOT the case that both a and b are finite. Here   
ln_ϕ′(y) = log(2) - 2 log(1 + y)                if -∞ = a < b < ∞
ln_ϕ′(y) = log(2) - 2 log(1 - y)                if -∞ < a < b = ∞
ln_ϕ′(y) = log(1 + y²) - 2 log(1 - y²)          if -∞ = a < b = ∞ 
"""
function ln_ϕ′_fn(T, y::AbstractFloat, a::AbstractFloat, b::AbstractFloat) 
    @assert a < b
    T_1 = convert(T, 1)
    T_2 = convert(T, 2)
    @assert -T_1 ≤ y ≤ T_1
    if -Inf < a && b < Inf
        throw(DomainError([a,b], "a and b are both finite"))
    elseif -Inf == a && b < Inf
        out = log(T_2) - T_2 * log(T_1 + y)
        return(out)
    elseif -Inf < a && b == Inf
        out = log(T_2) - T_2 * log(T_1 - y)
        return(out)
    else
        out = log(T_1 + y^2) - T_2 * log(T_1 - y^2)
        return(out)
    end
end


"""
# nodes_weights_support_fn & old_nodes_weights_support_fn
# produce the same results for the examples I have tried.
"""
function old_nodes_weights_support_fn(T, f_fn::Function, a, b, r::Integer)
    if T == Float64
        y, w = gausslegendre(r)
    else
        y, w = gauss(T, r)
    end
    nodes_support = ϕ_fn.(T, y, a, b) 
    weights_support = w .* f_fn.(nodes_support) .* ϕ′_fn.(T, y, a, b)
    [nodes_support, weights_support] 
end


"""
nodes_support, weights_support = nodes_weights_support_fn(T, f_fn, a, b, r)

Returns (x₁,...,xᵣ) and (w₁, ..., wᵣ), where 
    r                                         b
    ∑ wᵢ g(xᵢ) is the approximation used for  ∫ g(x) f(x) dx. 
   i=1                                        a             
Here f is the specified nonnegative weight function. 
This approximation is obtained by first approximating 
   1                          r
   ∫ h(y) dy       by         ∑ ξᵢ h(yᵢ),
  -1                         i=1 
where y₁,...,yᵣ are the Gauss Legendre nodes and ξ₁,...,ξᵣ 
are the corresponding weights. Since
   b                1
   ∫ g(x) f(x) dx = ∫ h(y) dy, where h(y) = g(ϕ(y)) f(ϕ(y)) ϕ′(y),        
   a               -1     
this leads us to approximate
   b                   r                                 r
   ∫ g(x) f(x) dx  by  ∑ ξᵢ g(ϕ(yᵢ)) f(ϕ(yᵢ)) ϕ′(yᵢ)  =  ∑ wᵢ g(xᵢ),             
   a                  i=1                               i=1                           i=1
where 
   wᵢ = ξᵢ ϕ′(yᵢ) f(ϕ(yᵢ)) and xᵢ = ϕ(yᵢ) (i=1,...,r).
"""
function nodes_weights_support_fn(T, f_fn::Function, a, b, r::Integer)
    if T == Float64
        y, w = gausslegendre(r)
    else
        y, w = gauss(T, r)
    end
    nodes_support = ϕ_fn.(T, y, a, b) 
    weights_support = w .* (f_fn.(nodes_support) .* ϕ′_fn.(T, y, a, b))
    [nodes_support, weights_support] 
end


"""
nodes_support, lnweights_support = nodes_lnweights_support_fn(T, f_fn, a, b, r)

This function is based on the function nodes_weights_support_fn,
which can result in NaN's, when using Double64 arithmetic,
for some of the weights w₁, ...,wᵣ. 
The present function solves this problem by computing
(log(w₁), ..., log(wᵣ)) instead.    
Returns (x₁,...,xᵣ) and (log(w₁), ..., log(wᵣ)), where 
    r                                         b
    ∑ wᵢ g(xᵢ) is the approximation used for  ∫ g(x) f(x) dx. 
   i=1                                        a             
Here f is the specified nonnegative weight function. 
This approximation is obtained by first approximating 
   1                          r
   ∫ h(y) dy       by         ∑ ξᵢ h(yᵢ),
  -1                         i=1 
where y₁,...,yᵣ are the Gauss Legendre nodes and ξ₁,...,ξᵣ 
are the corresponding weights. Since
   b                1
   ∫ g(x) f(x) dx = ∫ h(y) dy, where h(y) = g(ϕ(y)) f(ϕ(y)) ϕ′(y),        
   a               -1     
this leads us to approximate
   b                   r                                 r
   ∫ g(x) f(x) dx  by  ∑ ξᵢ g(ϕ(yᵢ)) f(ϕ(yᵢ)) ϕ′(yᵢ)  =  ∑ wᵢ g(xᵢ),             
   a                  i=1                               i=1                          
where 
   wᵢ = ξᵢ ϕ′(yᵢ) f(ϕ(yᵢ)) and xᵢ = ϕ(yᵢ) (i=1,...,r).
"""
function nodes_lnweights_support_fn(T, lnf_fn::Function, a, b, r::Integer)
    if T == Float64
        y, w = gausslegendre(r)
    else
        y, w = gauss(T, r)
    end
    nodes_support = ϕ_fn.(T, y, a, b) 
    lnweights_support = log.(w) + ln_ϕ′_fn.(T, y, a, b) + lnf_fn.(nodes_support) 
    [nodes_support, lnweights_support] 
end




"""
scaled_chi_pdf = scaled_chi_pdf_fn(T, x, m)

Returns the scaled chi probability density function (pdf) 
with m degrees of freedom, evaluated at x.
"""
function scaled_chi_pdf_fn(T, x::AbstractFloat, m::Integer)
    T_0 = convert(T,0)
    @assert m > T_0
    @assert x ≥ T_0
    T_2 = convert(T,2)
    tmp = m^(m/T_2) / gamma(m/T_2)
    term = tmp / T_2^((m/T_2) - 1)
    term * x^(m-1) * exp(- m * x^2 / T_2)
end




"""
ln_scaled_chi_pdf = ln_scaled_chi_pdf_fn(T, x, m)

Returns log of the scaled chi probability density function (pdf) 
with m degrees of freedom, evaluated at x.
"""
function ln_scaled_chi_pdf_fn(T, x::AbstractFloat, m::Integer)
    T_0 = convert(T,0)
    @assert m > T_0
    @assert x ≥ T_0
    T_1 = convert(T,1)
    T_2 = convert(T,2)
    T_m = convert(T,m)
    tmp1 = (T_m / T_2) * log(T_m)
    tmp2 = lgamma(T_m / T_2)
    tmp3 = ((T_m/T_2) - T_1) * log(T_2)
    tmp4 = (T_m - T_1) * log(x)
    tmp5 = T_m * x^2 / T_2
    tmp1 - tmp2 - tmp3 + tmp4 - tmp5
end



"""
lnf_chemistry = lnf_chemistry_fn(T, x)

Returns log of the weight function considered by Gaustschi (1983),
evaluated at x. This weight function is     
  f(x) = exp(-x³/3) for x > 0
       = 0  otherwise.
Reference
Gautschi, W. (1983) How and how not to check Gaussian quadrature
formulae. BIT, 23, 209-216.       
"""
function lnf_chemistry_fn(T::Type, x::AbstractFloat)
    T_0 = convert(T,0)
    @assert x ≥ T_0
    T_3 = convert(T,3)
    - convert(T, x)^T_3 / T_3
end



"""
lnf_hermite = lnf_hermite_fn(x)

Returns log of the Hermite weight function),
evaluated at x. This weight function is     
  f(x) = exp(-x²).     
"""
function lnf_hermite_fn(x::AbstractFloat)
    - x * x
end



"""
lnf_laguerre = lnf_laguerre_fn(x, α)

Returns log of the Generalized Laguerre weight function),
evaluated at x. This weight function is     
  f(x) = x^α exp(-x) for x > 0
       = 0  otherwise.    
"""
function lnf_laguerre_fn(x::AbstractFloat, α::AbstractFloat)
    @assert x ≥ convert(BigFloat, 0)
    @assert α > convert(BigFloat, -1)
    α * log(x) - x
end



"""
inner_prod_fns = inner_prod_fns_fn(u_fn, v_fn, nodes_support, weights_support)

Returns an approximation to the inner product of two functions
u: ℝ → ℝ and v: ℝ → ℝ. This inner product is defined to be
    ∞
    ∫ u(x) v(x) f(x) dx,                                  
   -∞
where f denotes the specified nonnegative weight function. 
We suppose that the support of the weight function f
is an interval with lower endpoint a and upper endpoint b, 
-∞ ≤ a < b ≤ ∞, so that this integral becomes
    b
    ∫ g(x) f(x) dx,  where the function g = uv.                                    
    a
"""
function inner_prod_fns_fn(u_fn, v_fn, nodes_support, weights_support)
    tmp1_vec = weights_support .* u_fn.(nodes_support) 
    tmp2_vec = v_fn.(nodes_support)
    dot(tmp1_vec, tmp2_vec)
end



"""
inner_prod_fns_nonan = inner_prod_fns_nonan_fn(u_fn, v_fn, nodes_support, lnweights_support)

This function is based on the function inner_prod_fns_fn,
which uses the weights w₁, ...,wᵣ. However, the computation of 
these weights using the function nodes_weights_support_fn 
can result in NaN's for some of weights, when using Double64 
arithmetic. The present function solves this problem by using 
(log(w₁), ..., log(wᵣ)) instead.    
Returns an approximation to the inner product of two functions
u: ℝ → ℝ and v: ℝ → ℝ. This inner product is defined to be
    ∞
    ∫ u(x) v(x) f(x) dx,                                  
   -∞
where f denotes the specified nonnegative weight function. 
We suppose that the support of the weight function f
is an interval with lower endpoint a and upper endpoint b, 
-∞ ≤ a < b ≤ ∞, so that this integral becomes
    b
    ∫ g(x) f(x) dx,  where the function g = uv.                                    
    a
The approximation to the inner product is computed so at to 
avoid resulting in a NaN.    
"""
function inner_prod_fns_nonan_fn(u_fn, v_fn, nodes_support, lnweights_support)
    u_vec = u_fn.(nodes_support)
    v_vec = v_fn.(nodes_support)
    tmp1_vec = sign.(u_vec) .* sign.(v_vec)
    tmp2_vec = exp.(lnweights_support .+ log.(abs.(u_vec)) .+ log.(abs.(v_vec)))
    dot(tmp1_vec, tmp2_vec)
end

   
