# This Julia script is part of the module 
# CustomGaussQuadrature.jl
# This script consists of
# code to compute the custom-made Gauss quadrature 
# nodes and weights using a moment-based method via 
# moment determinants to compute the coefficients
# in the three-term recurrence relation, using
# the type T, which is specified by the user. 
# This is followed by the use of the package
# GenericLinearAlgebra.jl


"""
moment = moment_fn(T, which_f, r)

This function computes moment of order r
for the weight function specified by which_f
The r'th moment is
      ∞
 μᵣ = ∫ xʳ f(x) dx,
     −∞
for any nonnegative integer r. The standard notation
for this moment is μᵣ′, so the notation μᵣ is not standard.
However, μᵣ is a convenient notation and accords with that
used by 
Gautschi, W. (2004) Orthogonal Polynomials, 
Computation and Approximation. Oxford University Press, Oxford.
These are the moments needed to compute
the custom-made Gauss quadrature rule for
n nodes.
The input which_f has the following 3 components:
(i)  name
(ii) support specified by a 2-vector of the endpoints 
of the interval, with finite endpoints specified by 
integers (which are later converted to the appropriate 
floating-point type)
(iii) parameter vector (if any)
This function accepts the following what_if's as inputs:
which_f = ["scaled chi pdf", [0,Inf], m]
which_f = ["Hermite", [-Inf, Inf]]
which_f = ["Generalized Laguerre", [0, Inf], α_GGL]     
which_f = ["chemistry example", [0, Inf]]
which_f = ["Legendre", [-1, 1]]
"""
function moment_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
  @assert r ≥ 0
  T_2 = convert(T, 2)

  if which_f[1] == "scaled chi pdf"
    if r == 0
        return(convert(T, 1))
    end
    m = which_f[3]
    T_m = convert(T, m)
    T_r = convert(T, r)
    term1 = (T_r/T_2) * log(T_2 / T_m)
    # 3 October 2024: loggamma was not recognised
    # but lgamma was. So I used
    # term2 = lgamma(((T_r + T_m)/T_2)))
    # term3 = lgamma(T_m/T_2)
    # However, when testing the package 
    # CustomGaussQuadrature on 7 February 2025, 
    # there was a warning that lgammm(x::Real) is
    # deprecated and should be replaced
    # by (logabsgamma(x))[1]
    term2 = (logabsgamma((T_r + T_m)/T_2))[1]
    term3 = (logabsgamma(T_m/T_2))[1]
    moment = exp(term1 + term2 - term3)
    return(moment)

  elseif which_f[1] == "Hermite"
    if r == 0
        return(sqrt(convert(T, π)))
    end
    if isodd(r)
        moment = convert(T, 0)
    else
        num  = convert(T, r + 1)
        moment = gamma(num / T_2)
    end 
    return(moment)

  elseif which_f[1] == "Generalized Laguerre"
    α_GGL = which_f[3]
    term = convert(T, r + α_GGL + 1)
    moment = gamma(term)
    return(moment)
    
  elseif which_f[1] == "Legendre"
    if isodd(r)
        moment = convert(T, 0)
    else
        denom = convert(T, r + 1)
        moment = T_2 / denom
    end
    return(moment)

  elseif which_f[1] ==  "chemistry example"
    T_1 = convert(T, 1)
    T_3 = convert(T, 3)
    T_r = convert(T, r)
    term1 = T_3^((T_r - T_2) / T_3)
    term2 = gamma((T_r + T_1) / T_3)
    moment = term1 * term2
    return(moment)

  else
    # DomainError means that the argument to a function
    # or constructor does not lie in the valid domain
    throw(DomainError(which_f, "invalid argument"))
  end
end



"""
μ_offsetvec = μ_offsetvec_fn(T, which_f, n)

This function places the moments of orders
0, ... , 2n - 1
for the weight function specified by which_f
in an offset vector, whose starting index is 0. 
The r'th moment is
      ∞
 μᵣ = ∫ xʳ f(x) dx,
     −∞
for any nonnegative integer r. The standard notation
for this moment is μᵣ′, so the notation μᵣ is not standard.
However, μᵣ is a convenient notation and accords with that
used by 
Gautschi, W. (2004) Orthogonal Polynomials, 
Computation and Approximation. Oxford University Press, Oxford.
These are the moments needed to compute
the custom-made Gauss quadrature rule for
n nodes.
The input which_f has the following 3 components:
(i)  name
(ii) support specified by a 2-vector of the endpoints 
of the interval, with finite endpoints specified by 
integers (which are later converted to the appropriate 
floating-point type)
(iii) parameter vector (if any)
"""
function μ_offsetvec_fn(::Type{T}, which_f, n::Integer) where {T<:AbstractFloat}
  @assert n ≥ 1
  μ_offsetvec = OffsetVector(zeros(T,2*n), 0:(2*n - 1))
  for r in 0:(2*n-1)
    μ_offsetvec[r] = moment_fn(T, which_f, r)
  end
  μ_offsetvec
end



"""
Δ = Δ_fn(μ_offsetvec, k, n)

This function computes the Hankel determinant
of order k, where 0 ≤ k ≤ n, as defined by (2.1.1) 
on page 53 of
Gautschi, W. (2004) Orthogonal Polynomials:
Computation and Approximation.
Here n is the number of Gauss quadrature nodes.
"""  
function Δ_fn(μ_offsetvec, k::Integer, n::Integer) 
 @assert n ≥ 2 
 @assert length(μ_offsetvec) == 2*n
 @assert 0 ≤ k ≤ n  
 T = typeof(μ_offsetvec[0])  
 if k == 0
  return(convert(T,1))
 end
 Hankel_mat = zeros(T, k, k)
  for i in 1:k
    for j in 1:k
      r = (j - 1) + (i - 1)
      Hankel_mat[i, j] = μ_offsetvec[r]
    end
  end
  det(Symmetric(Hankel_mat))
end



"""
Δ′ = Δ′_fn(μ_offsetvec, k, n)

This function computes the determinant of order k, 
where 0 ≤ k ≤ n, as defined by (2.1.2) on page 53 of
Gautschi, W. (2004) Orthogonal Polynomials:
Computation and Approximation.
Here n is the number of Gauss quadrature nodes.
"""  
function Δ′_fn(μ_offsetvec, k::Integer, n::Integer) 
  @assert n ≥ 2 
  @assert length(μ_offsetvec) == 2*n  
  @assert 0 ≤ k ≤ n  
  T = typeof(μ_offsetvec[0]) 

  if k == 0
    return(convert(T, 0))
  elseif k == 1
    return(μ_offsetvec[1])
  end

  Hankel_prime_part1_mat = zeros(T, k, (k-1))
  for i in 1:k
    for j in 1:(k-1)
      r = (j - 1) + (i - 1)
      Hankel_prime_part1_mat[i, j] = μ_offsetvec[r]
    end
  end
  
  Hankel_prime_last_col = zeros(T, k)
  for i in 1:k
    r = i + k - 1
    Hankel_prime_last_col[i] = μ_offsetvec[r]
  end
  
  Hankel_prime_mat = hcat(Hankel_prime_part1_mat, Hankel_prime_last_col)

  det(Hankel_prime_mat)
end



"""
Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n)

This function computes the vector [Δ₀, ..., Δₙ],
where n is the number of Gauss quadrature nodes.
This is an offset vector with the index to its successive
elements starting at 0.
 """ 
function Δ_offsetvec_fn(μ_offsetvec, n::Integer)
  @assert n ≥ 2
  @assert length(μ_offsetvec) == 2*n
  T = typeof(μ_offsetvec[0])
  out_offsetvec = OffsetVector(zeros(T, (n+1)), 0:n)
  for k in 0:n
      out_offsetvec[k] = Δ_fn(μ_offsetvec, k, n)
  end
  out_offsetvec
end


"""
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n)

This function computes the vector [Δ₀′, ..., Δₙ′],
where n is the number of Gauss quadrature nodes.
This is an offset vector with the index to its successive
elements starting at 0.
 """ 
function Δ′_offsetvec_fn(μ_offsetvec, n::Integer)
  @assert n ≥ 2
  @assert length(μ_offsetvec) == 2*n
  T = typeof(μ_offsetvec[0])
  out_offsetvec = OffsetVector(zeros(T, (n+1)), 0:n)
  for k in 0:n
      out_offsetvec[k] = Δ′_fn(μ_offsetvec, k, n)
  end
  out_offsetvec
end


"""
α_offsetvec = α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n)

This function computes the vector [α₀, α₁, ..., αₙ₋₁],
where n is the number of Gauss quadrature nodes.
This is an offset vector with the index to its successive
elements starting at 0.
αₖ is defined by (2.1.4) on p.54 of 
Gautschi, W. (2004) Orthogonal Polynomials: 
Computation and Approximation.
 """ 
function α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n::Integer)
  @assert n ≥ 2
  @assert length(Δ_offsetvec) == n+1
  @assert length(Δ′_offsetvec) == n+1
  T = typeof(Δ_offsetvec[0])
  out_offsetvec = OffsetVector(zeros(T, n), 0:(n-1))
  for k in 0:(n-1)
      term1 = Δ′_offsetvec[k+1] / Δ_offsetvec[k+1]
      term2 = Δ′_offsetvec[k] / Δ_offsetvec[k]
      out_offsetvec[k] = term1 - term2
  end
  out_offsetvec
end



"""
β_vec = β_vec_fn(Δ_offsetvec, n)

This function computes the vector [β₁, ..., βₙ],
where n is the number of Gauss quadrature nodes.
βₖ is defined by by (2.1.5) on p.54 of
Gautschi, W. (2004) Orthogonal Polynomials:
Computation and Approximation.
 """ 
function β_vec_fn(Δ_offsetvec, n::Integer)
  @assert n ≥ 2
  @assert length(Δ_offsetvec) == n+1
  T = typeof(Δ_offsetvec[0])
  out_vec =zeros(T, (n-1))
  for k in 1:(n-1)
    out_vec[k] = Δ_offsetvec[k+1] * Δ_offsetvec[k-1] / Δ_offsetvec[k]^2
  end
  out_vec
end



# 2 Oct 24: I should be able to delete the following
# function at some stage.
function α_offsetvec_β_vec_fn(T, which_f, n::Integer)
  μ_offsetvec = μ_offsetvec_fn(T, which_f, n) 
  Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n)
  Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n)
  α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n)
  β_vec = β_vec_fn(μ_offsetvec, Δ_offsetvec, n)
  [α_offsetvec, β_vec]
end

"""
a_vec, b_vec, μ₀ = a_vec_b_vec_μ₀_fn(T, which_f, n)

This function computes a_vec, b_vec and μ₀ using the
type of arithmetic T (usually a BigFloat with a globally
specified precision), which_f which has the following 
3 components:
(i)  name
(ii) support specified by a 2-vector of the endpoints 
of the interval, with finite endpoints specified by 
integers (which are later converted to the appropriate 
floating-point type)
(iii) parameter vector (if any)
and the number n of Gauss quadrature nodes.
"""
function a_vec_b_vec_μ₀_fn(T, which_f, n::Integer)
  @assert n ≥ 2
  μ_offsetvec = μ_offsetvec_fn(T, which_f, n) 
  μ₀ =  μ_offsetvec[0]
  Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n)
  β_vec = β_vec_fn(Δ_offsetvec, n)
  T_0 = convert(T, 0)
  @assert all(x -> x >= T_0, β_vec)
  b_vec = sqrt.(β_vec)
  Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n)
  α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n)
  # The parent function converts an Offset Array 
  # into an ordinary Array
  a_vec = parent(α_offsetvec)
  [a_vec, b_vec, μ₀]
end



"""
a_vec, b_vec, μ₀, nbits = step1a_fn(which_f, n)

This function is the first part of Step 1.
It chooses the initial number of bits
of precision for BigFloat such that all of the 
elements of a_vec and b_vec can be computed.
Of course, this does not mean that a_vec and 
b_vec are computed to sufficient precision for
the computation of the Gauss quadrature 
coefficients to sufficient precision.  
"""
function step1a_fn(which_f, n::Integer)
  @assert n ≥ 2
  a_vec, b_vec, μ₀ =
  try 
      setprecision(BigFloat, 80, base=2)
      a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
  catch 
      try
          setprecision(BigFloat, 106, base=2)
          a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
      catch
          try
              setprecision(BigFloat, 132, base=2)
              a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
          catch
              try
                  setprecision(BigFloat, 158, base=2)
                  a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
              catch
                  try
                      setprecision(BigFloat, 184, base=2)
                      a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                  catch
                      try
                          setprecision(BigFloat, 210, base=2)
                          a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                      catch
                          try
                              setprecision(BigFloat, 236, base=2)
                              a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                          catch
                              try
                                  setprecision(BigFloat, 262, base=2)
                                  a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                              catch
                                  try
                                      setprecision(BigFloat, 288, base=2)
                                      a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                                  catch
                                      try
                                          setprecision(BigFloat, 340, base=2)
                                          a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                                      catch
                                          try
                                              setprecision(BigFloat, 366, base=2)
                                              a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                                          catch
                                              try
                                              setprecision(BigFloat, 392, base=2)
                                              a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
                                              finally
                                              throw(DomainError(n, " too large"))
                                              end
                                          end
                                      end
                                  end
                              end
                          end 
                      end
                  end
              end
          end
      end
  end
nbits = precision(BigFloat)

# The precision of BigFloat is set globally.
# Consequently, this precision is reset to
# its default value before exitting this
# function. 
setprecision(BigFloat, 256, base=2)

[a_vec, b_vec, μ₀, nbits]
end



"""
a_vec, b_vec, μ₀, nbits = step1_fn(which_f, n)

This function carries our Step 1 which is
the computation of a_vec, b_vec and μ₀
to sufficient BigFloat precision so that 
the computation of the Gauss quadrature nodes
and weights can be carried out in Step 2 to
full Float64 precision. This sufficient 
BigFloat precision is given by nbits. 
"""
function step1_fn(which_f, n::Integer, epsilon=1.0e-18)
  @assert n ≥ 2

  a_vec, b_vec, μ₀, nbits = step1a_fn(which_f, n)
    
  nbits_new = nbits+26
  setprecision(BigFloat, nbits_new, base=2)
  a_vec_new, b_vec_new, μ₀_new = a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
 
  abs_error_a = convert(Vector{Float64}, abs.(a_vec - a_vec_new))
  max_abs_error_a = maximum(abs_error_a)
  rel_error_b = convert(Vector{Float64}, (abs.((b_vec - b_vec_new) ./ b_vec_new)))
  max_rel_error_b = maximum(rel_error_b)
    
  while (max_abs_error_a > epsilon || max_rel_error_b > epsilon) && nbits_new <= 444
    nbits = nbits_new
    a_vec, b_vec, μ₀ = a_vec_new, b_vec_new, μ₀_new
    nbits_new = nbits + 26
    setprecision(BigFloat, nbits_new, base=2)
    a_vec_new, b_vec_new, μ₀_new = a_vec_b_vec_μ₀_fn(BigFloat, which_f, n) 
    abs_error_a = convert(Vector{Float64}, abs.(a_vec - a_vec_new))
    max_abs_error_a = maximum(abs_error_a)
    rel_error_b = convert(Vector{Float64}, (abs.((b_vec - b_vec_new) ./ b_vec_new)))
    max_rel_error_b = maximum(rel_error_b)
  end

  # The precision of BigFloat is set globally.
  # Consequently, this precision is reset to
  # its default value before exitting this
  # function. 
  setprecision(BigFloat, 256, base=2)

  [a_vec_new, b_vec_new, μ₀_new, nbits_new]
end



"""
nodes, weights = step2_fn(T, which_f, n, a_vec, b_vec, μ₀)

This function carries our Step 2 which is the computation of the
nodes and weights from the eigenvalues and eigenvectors
of the symmetric tridiagonal Jacobi matrix Jₙ with main diagonal
a_vec and subdiagonal b_vect. 
This computation is carried out using the type T of floating-point 
arithmetic. 
For flexibility, the resulting nodes and weights are outputted
as vectors whose elements are of the same type. 
Usually the type T is Double64 or BigFloat with a globally set
precision. 
"""
function step2_fn(T, which_f, n::Integer, a_vec, b_vec, μ₀)
  @assert n ≥ 2
  a_vec_converted = convert(Vector{T}, a_vec)
  b_vec_converted = convert(Vector{T}, b_vec)
  μ₀_converted = convert(T, μ₀)
  Jₙ = SymTridiagonal(a_vec_converted, b_vec_converted)
  eigenvalues, eigenvectors = eigen(Jₙ)
  nodes = eigenvalues
  endpts = which_f[2]
  @assert convert(T, endpts[1]) ≤ nodes[1]
  @assert nodes[n] ≤ convert(T, endpts[2])
  weights = zeros(T,n)
  for j in 1:n
    weights[j] = eigenvectors[1, j]^2 / dot(eigenvectors[:,j], eigenvectors[:,j])
  end
  weights = μ₀_converted * weights
  [nodes, weights]
end



"""
nodes, weights = custom_gauss_quad_fn(T, which_f, n)

This function carries out Step 1, followed by Step 2, which is 
carried out using the type T of floating-point arithmetic. 
For flexibility, the resulting nodes and weights are outputted
as vectors whose elements are of the same type. 
Usually the type T is Double64 or BigFloat with a globally set
precision. 
"""
function custom_gauss_quad_fn(T, which_f, n::Integer)
  @assert n ≥ 1
  if n ==1
    μ₀, μ₁ = μ_offsetvec_fn(T, which_f, n)
    nodes = μ₁ / μ₀
    weights = μ₀
    return([nodes, weights])
  end
  a_vec, b_vec, μ₀ = step1_fn(which_f, n)
  nodes, weights = step2_fn(T, which_f, n, a_vec, b_vec, μ₀)
  [nodes, weights]
end



"""
nodes, weights = custom_gauss_quad_all_fn(which_f, n, upto_n=false, extra_check=false)

This function computes the custom-made Gauss quadrature nodes and 
weights, for n nodes and the weight function specified by the input 
which_f with the following 3 components:
(i)  name
(ii) support specified by a 2-vector of the endpoints 
of the interval, with finite endpoints specified by 
integers (which are later converted to the appropriate 
floating-point type)
(iii) parameter vector (if any)
This function accepts the following what_if's as inputs:
which_f = ["scaled chi pdf", [0,Inf], m]
which_f = ["Hermite", [-Inf, Inf]]
which_f = ["Generalized Laguerre", [0, Inf], α_GGL]     
which_f = ["chemistry example", [0, Inf]]
which_f = ["Legendre", [-1, 1]]
The Boolean variable upto_n takes the value true if the Gauss quadrature 
rules are computed for number of nodes q from 1 up and including n.
By default, this variable is false, so that
only the Gauss quadrature rule for number of nodes n is computed.
The Boolean variable extra_check takes the value true if 26 bits of 
precision are added to the number of bits of precision (a) used in the 
final stage of Step 1 and (b) Step 2, followed by comparison of this 
more precise computation of the nodes and weights with the less 
precise computation of the nodes and weights.
By default, this variable is false, so that this extra check on the
precision of the nodes and weights is not carried out. 
"""
function custom_gauss_quad_all_fn(which_f, n::Integer, upto_n::Bool=false, extra_check::Bool=false)
  @assert n ≥ 1
  if n == 1
    T = Double64  
    μ₀, μ₁ = μ_offsetvec_fn(T, which_f, n)
    nodes = μ₁ / μ₀
    weights = μ₀
    return([nodes, weights])
  end
  a_vec, b_vec, μ₀, nbits = step1_fn(which_f, n)
  nodes, weights = step2_fn(Double64, which_f, n, a_vec, b_vec, μ₀)
  if extra_check == true
    setprecision(BigFloat, (nbits+26), base=2)
    a_vec_better, b_vec_better, μ₀_better = step1_fn(which_f, n)
    setprecision(BigFloat, 132, base=2)
    nodes_better, weights_better = step2_fn(BigFloat, which_f, n, a_vec_better, b_vec_better, μ₀_better)
    nodes = convert(Vector{BigFloat}, nodes)
    weights = convert(Vector{BigFloat}, weights)
    abs_error_nodes = convert(Vector{Float64}, abs.(nodes - nodes_better));
    max_abs_error_nodes = maximum(abs_error_nodes)
    abs_rel_error_weights = convert(Vector{Float64}, (abs.((weights - weights_better) ./ weights_better)))
    max_abs_rel_error_weights = maximum(abs_rel_error_weights)
    a_vec = a_vec_better
    b_vec = b_vec_better
    nodes = nodes_better
    weights = weights_better
  end

  if upto_n == false
    if extra_check == true
      # The precision of BigFloat is set globally.
      # Consequently, this precision is reset to
      # its default value before exitting this
      # function. 
      setprecision(BigFloat, 256, base=2)
      return([nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights])
    end
    # The precision of BigFloat is set globally.
    # Consequently, this precision is reset to
    # its default value before exitting this
    # function. 
    setprecision(BigFloat, 256, base=2)
    return([nodes, weights])
  end
  nodes_upto_n = Any[]
  weights_upto_n = Any[]
  T = Double64
  μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1)
  nodes_1 = μ₁ / μ₀
  weights_1 = μ₀
  push!(nodes_upto_n, nodes_1)
  push!(weights_upto_n, weights_1)
  for q in 2:(n-1)
    a_vec_q = a_vec[1:q]
    b_vec_q = b_vec[1:q]
    nodes_q, weights_q = step2_fn(Double64, which_f, q, a_vec_q, b_vec_q, μ₀)
    push!(nodes_upto_n, nodes_q)
    push!(weights_upto_n, weights_q)
  end
  push!(nodes_upto_n, nodes)
  push!(weights_upto_n, weights)
  if extra_check == true
    # The precision of BigFloat is set globally.
    # Consequently, this precision is reset to
    # its default value before exitting this
    # function. 
    setprecision(BigFloat, 256, base=2)
    return([nodes_upto_n, weights_upto_n, max_abs_error_nodes, max_abs_rel_error_weights])
  end

  # The precision of BigFloat is set globally.
  # Consequently, this precision is reset to
  # its default value before exitting this
  # function. 
  setprecision(BigFloat, 256, base=2)

  [nodes_upto_n, weights_upto_n]
end







