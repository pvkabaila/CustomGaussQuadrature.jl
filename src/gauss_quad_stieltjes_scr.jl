# This Julia script is part of the module 
# CustomGaussQuadrature.jl
#
# This script consists of Julia functions to compute 
# the Gauss quadrature rule with n nodes using the 
# following two steps:
#
# Step 1: Compute the recursion coefficients 
# [α₀, α₁, ..., αₙ₋₁] and [β₁, ..., βₙ] of the 
# three-term recurrence relation using the 
# Stieltjes procedure, where the inner product 
# of two functions is found using the script 
# inner_prod_fns_scr.jl. The aim of 
# Step 1 is to compute approximations to the 
# vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and 
# [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁] with maximum
# absolute errors and maximum relative errors,
# respectively, bounded above by 10⁻¹⁸.
#
# Step 2: Compute the Gauss quadrature rule with n nodes
# from the eigenvalues and eigenvectors of the symmetric 
# tridiagonal Jacobi matrix Jₙ with main diagonal
# [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and subdiagonal 
# [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], using  
# Double64 arithmetic. The computation of these 
# eigenvalues and eigenvectors is carried out using
# the package GenericLinearAlgebra.jl

"""
a_vec, b_vec = stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_support, lnweights_support)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and 
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], where n is the number 
of Gauss quadrature nodes (n ≥ 2). αₖ and βₖ are defined by 
  αₖ = (x πₖ, πₖ)ₐ / (πₖ, πₖ)ₐ and  βₖ = (πₖ, πₖ)ₐ / (πₖ₋₁, πₖ₋₁)ₐ 
where (u,v)ₐ denotes the discrete approximation               
  r                               b      
  ∑ wᵢ u(xᵢ) v(xᵢ)   to   (u,v) = ∫ u(x) v(x) f(x) dx.           
 i=1                              a                        
Here 
[x₁, ... , xᵣ] is the input nodes_support
[log(w₁), ... , log(wᵣ)] is the input lnweights_support
The type of arithmetic used in all computations is 
specified by the input T. If T is not appropriately 
chosen then the the vectors [a₁, ..., aₙ] and 
[b₁, ..., bₙ₋₁] may be insufficiently accurate and
may even include elements that are NaN's.
"""
function stieltjes_a_vec_b_vec_nonan_fn(T::Type, n::Integer, μ₀, nodes_support, lnweights_support)
    @assert n ≥ 2
    α_offsetvec = OffsetVector(zeros(T, n), 0:(n-1))
    β_vec = zeros(T, (n-1))
    sqnorms_offsetvec = OffsetVector(zeros(T, n), 0:(n-1))

    sqnorms_offsetvec[0] = μ₀
    poly_1 = Polynomial([convert(T,1)])
    poly_x = Polynomial([convert(T,0), convert(T,1)])
    num = inner_prod_fns_nonan_fn(poly_x, poly_1, nodes_support, lnweights_support)
    denom = sqnorms_offsetvec[0]
    α_offsetvec[0] = num / denom
    array_monomials = Any[]
    monomial = Polynomial([-α_offsetvec[0], convert(T,1)])
    push!(array_monomials, monomial)

    sqnorms_offsetvec[1] = inner_prod_fns_nonan_fn(array_monomials[1], array_monomials[1], 
    nodes_support, lnweights_support)
    β_vec[1] = sqnorms_offsetvec[1] / sqnorms_offsetvec[0]
    term = poly_x * array_monomials[1]
    num = inner_prod_fns_nonan_fn(term, array_monomials[1], 
        nodes_support, lnweights_support)
    denom = sqnorms_offsetvec[1]
    α_offsetvec[1] = num / denom    
    monomial = Polynomial([-α_offsetvec[1], convert(T,1)]) * array_monomials[1] - β_vec[1]
    push!(array_monomials, monomial)

    if n == 2
        a_vec = parent(α_offsetvec)
        b_vec = sqrt.(β_vec[1:(n-1)])
        return([a_vec, b_vec])
    end

    for k in 2:(n-1)
        sqnorms_offsetvec[k] = inner_prod_fns_nonan_fn(array_monomials[k], array_monomials[k], 
        nodes_support, lnweights_support)
        β_vec[k] = sqnorms_offsetvec[k] / sqnorms_offsetvec[(k-1)]
        term = poly_x * array_monomials[k]
        num = inner_prod_fns_nonan_fn(term, array_monomials[k], 
            nodes_support, lnweights_support)
        denom = sqnorms_offsetvec[k]
        α_offsetvec[k] = num / denom    
        monomial = Polynomial([-α_offsetvec[k], convert(T,1)]) * array_monomials[k] - β_vec[k] * array_monomials[(k-1)]
        push!(array_monomials, monomial)
    end
   
    a_vec = parent(α_offsetvec)
    b_vec = sqrt.(β_vec[1:(n-1)])

    [a_vec, b_vec]
end


# Materialize (convert to a concrete value of type T) a Stieltjes scalar input,
# either from a user-supplied specification or from a callback.
function stieltjes_materialize_typed_scalar_fn(value_or_fn, T::Type, which_f)
  if value_or_fn isa Function
    try
      return value_or_fn(T, which_f)
    catch err
      if err isa MethodError
        return value_or_fn(T)
      end
      rethrow()
    end
  end
  return materialize_scalar_spec_fn(T, value_or_fn)
end


function stieltjes_a_vec_b_vec_choosenbits_core_fn(n, make_μ₀_fn::Function, make_lnf_fn::Function, a, b, r)
  setprecision(BigFloat, 256, base=2)
  T = BigFloat
  nbits = 256
  # Start from a 256-bit BigFloat run and treat it as the reference solution.
  # Later branches only decide whether this precision is already sufficient or
  # whether we must move up to 512 bits.
  # Build lnf_fn only after choosing T, so every downstream computation
  # sees a closure whose constants and branches are created for that T.
  # The support endpoints a and b are still passed in as raw values; the
  # support-quadrature helper converts them to this local T, so endpoint
  # storage does not bias the BigFloat-versus-Double64 comparison.
  μ₀ = make_μ₀_fn(T)
  lnf_fn = make_lnf_fn(T)
  nodes_256, lnweights_256 = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_256, b_vec_256 = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_256, lnweights_256)

  epsilon1 = 1.0e-22

  T = Double64
  # Cheap consistency check: if Double64 reproduces the 256-bit coefficients
  # to within the target tolerances and without NaNs, then 256 bits is already
  # more than enough and we can return the 256-bit result.
  μ₀ = make_μ₀_fn(T)
  lnf_fn = make_lnf_fn(T)
  nodes_106, lnweights_106 = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_106, b_vec_106 = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_106, lnweights_106)
  if any(isnan, a_vec_106)==false && any(isnan, b_vec_106)==false
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_106 - a_vec_256))
    max_abs_error_a = maximum(abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_106 - b_vec_256) ./ b_vec_256)))
    max_rel_error_b = maximum(rel_error_b)
    if max_abs_error_a < epsilon1 && max_rel_error_b < epsilon1
      return([a_vec_256, b_vec_256, nbits])
    end
  end

  setprecision(BigFloat, 224, base=2)
  T = BigFloat;
  # Second check in the opposite direction: compare a slightly lower BigFloat
  # precision against the 256-bit reference. Agreement here indicates that the
  # 256-bit result is internally stable, so we keep the 256-bit output.
  μ₀ = make_μ₀_fn(T)
  lnf_fn = make_lnf_fn(T)
  nodes_224, lnweights_224 = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_224, b_vec_224 = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_224, lnweights_224) 
  abs_error_a = convert(Vector{Float64}, abs.(a_vec_224 - a_vec_256))
  max_abs_error_a = maximum(abs_error_a)
  rel_error_b = convert(Vector{Float64}, 
  (abs.((b_vec_224 - b_vec_256) ./ b_vec_256)))
  max_rel_error_b = maximum(rel_error_b)
  if max_abs_error_a < epsilon1 && max_rel_error_b < epsilon1
    return([a_vec_256, b_vec_256, nbits])
  end
            
  setprecision(BigFloat, 512, base=2)
  T = BigFloat
  nbits = 512
  # If the 256-bit result fails both cheaper stability checks, recompute at
  # 512 bits and treat that higher-precision run as the new candidate output.
  μ₀ = make_μ₀_fn(T)
  lnf_fn = make_lnf_fn(T)
  nodes_512 , lnweights_512  = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_512 , b_vec_512  = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_512 , lnweights_512 ) 
  abs_error_a = convert(Vector{Float64}, abs.(a_vec_512 - a_vec_256))
  max_abs_error_a = maximum(abs_error_a)
  rel_error_b = convert(Vector{Float64}, 
  (abs.((b_vec_512 - b_vec_256) ./ b_vec_512)))
  max_rel_error_b = maximum(rel_error_b)
  if max_abs_error_a < epsilon1 && max_rel_error_b < epsilon1
    # Even though we evaluated 512 bits, the two runs are effectively the same;
    # return the 512-bit arrays because this branch has selected 512 bits as the
    # working precision for the final output.
    return([a_vec_512, b_vec_512, nbits])
  end
    
  setprecision(BigFloat, 480, base=2)
  T = BigFloat;
  # Final safeguard: 480-bit and 512-bit results should agree if 512 bits is
  # genuinely stable. A mismatch means the requested n is too large for this
  # precision ladder, so the function rejects the problem instead of returning
  # coefficients whose accuracy is uncertain.
  μ₀ = make_μ₀_fn(T)
  lnf_fn = make_lnf_fn(T)
  nodes_480, lnweights_480 = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_480, b_vec_480 = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_480, lnweights_480) 
  abs_error_a = convert(Vector{Float64}, abs.(a_vec_480 - a_vec_512))
  max_abs_error_a = maximum(abs_error_a)
  rel_error_b = convert(Vector{Float64}, 
  (abs.((b_vec_480 - b_vec_512) ./ b_vec_512)))
  max_rel_error_b = maximum(rel_error_b)
  if max_abs_error_a > epsilon1 || max_rel_error_b > epsilon1
    throw(DomainError(n, " too large"))  
  end    
  setprecision(BigFloat, 256, base=2)

  # Reaching this point means 512 bits passed the final stability check, so the
  # function returns the 512-bit recurrence coefficients together with nbits=512.
  [a_vec_512, b_vec_512, nbits]
end


"""
a_vec, b_vec, nbits = stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], where n is the number
of Gauss quadrature nodes (n ≥ 2).

The input lnf_typed_fn is expected to have the form
lnf_typed_fn(T, which_f, x). The μ₀ input may be either a numeric
value or a typed factory such as μ₀(T, which_f) or μ₀(T). The
algorithm first chooses the arithmetic type T and then constructs
the one-argument closure lnf_fn(x) = lnf_typed_fn(T, which_f, x)
that is used in the support quadrature and Stieltjes recurrence
computations. The support endpoints a and b are kept as raw values
until the support quadrature converts them to this same local T,
so their stored type does not influence which arithmetic branch is
used. The driver does not inspect which_f[3], so for user-defined
weights that third component may be a container holding several
parameter specifications.
"""
function stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)
  #
  # Read T -> (x -> lnf_typed_fn(T, which_f, x)) from left to right:
  # 1. The outer function takes one input, called T.
  # 2. After T has been chosen, it returns a new function of x.
  # 3. That returned function remembers the particular values of T and
  #    which_f that were in scope when it was created.
  #
  # This remembered behaviour is a closure. In other words, the inner
  # function x -> ... has access to T and which_f even though x is its
  # only explicit argument. So if the driver first chooses T = BigFloat,
  # the closure behaves like x -> lnf_typed_fn(BigFloat, which_f, x).
  # If the driver later chooses T = Double64, a different closure is
  # created, namely x -> lnf_typed_fn(Double64, which_f, x).
  #
  # The key point is scope: T is not looked up later in some ambient
  # global scope. It is fixed by the outer function call that creates the
  # closure, and the resulting one-argument lnf_fn is then passed into the
  # existing support-quadrature and Stieltjes code.
  stieltjes_a_vec_b_vec_choosenbits_core_fn(
    n,
    T -> stieltjes_materialize_typed_scalar_fn(μ₀, T, which_f),
    T -> (x -> lnf_typed_fn(T, which_f, x)),
    a,
    b,
    r,
  )
end


"""
a_vec, b_vec, nbits, r = stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=7, j_max=40, epsilon=1.0e-15)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], together with the chosen
BigFloat precision nbits and the final value of r.

The input lnf_typed_fn is expected to have the form
lnf_typed_fn(T, which_f, x). For each trial value of r, the
driver chooses the arithmetic type T and then creates the
one-argument closure lnf_fn(x) = lnf_typed_fn(T, which_f, x)
before carrying out the support quadrature and Stieltjes
computations. For user-defined weights, which_f[3] may therefore be
empty, scalar, or a container of several parameter specifications,
provided lnf_typed_fn and μ₀ know how to interpret it.
"""
function stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=7, j_max=80, epsilon=1.0e-15)
  @assert n ≥ 2
  # Start with two consecutive trial values of r so the driver can measure how
  # much the recurrence coefficients change when the support quadrature is refined.
  # Here j is the trial index, and each candidate r is generated from it by
  # r = j * (offset + n). Increasing j by 1 therefore advances to the next
  # larger candidate r in this search.
  j = 3
  r = j * (offset + n)
  a_vec, b_vec, nbits = 
    stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

  j = j + 1
  r = j * (offset + n)
  a_vec_new, b_vec_new, nbits_new = 
    stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

  # Compare successive runs in absolute error for a_vec and relative error for
  # b_vec. The output is accepted only when increasing r no longer changes the
  # coefficients by more than the target tolerance.
  abs_error_a = convert(Vector{Float64}, abs.(a_vec_new - a_vec))
  max_abs_error_a = maximum(abs_error_a)
  rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_new - b_vec) ./ b_vec_new)))
  max_rel_error_b = maximum(rel_error_b)

  # Keep increasing r until two consecutive coefficient sets agree to within
  # epsilon, or until j reaches the user-supplied safety limit j_max.
  # Equivalently, the loop advances through j = 3, 4, 5, ... and tests the
  # corresponding sequence of candidate r values defined by that formula.
  while (max_abs_error_a > epsilon || max_rel_error_b > epsilon)
    # Promote the most recent run to the reference solution before computing the
    # next candidate at a larger value of r.
    a_vec, b_vec, nbits = a_vec_new, b_vec_new, nbits_new
    if j == j_max
      throw(DomainError(j, " needed value of the r search index exceeds j_max"))
    end
    j = j + 1
    r = j * (offset + n)
    a_vec_new, b_vec_new, nbits_new = 
    stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)
    # The new r is therefore not chosen independently; it is exactly the value
    # associated with the updated integer search index j.

    # Re-evaluate the stopping test against the updated reference and candidate
    # pair; termination means the latest candidate is stable with respect to r.
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_new - a_vec))
    max_abs_error_a = maximum(abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_new - b_vec) ./ b_vec_new)))
    max_rel_error_b = maximum(rel_error_b)

  end

  setprecision(BigFloat, 256, base=2)
  # Return the last candidate, which is the first one whose coefficients agree
  # with the previous run to within the prescribed tolerance.
  [a_vec_new, b_vec_new, nbits_new, r]
end


"""
nodes, weights = stieltjes_step2_fn(n, μ₀, a_vec, b_vec, a, b)

This function carries our Step 2, which is the computation of 
the Gauss quadrature nodes and weights from the eigenvalues 
and eigenvectors of the symmetric tridiagonal Jacobi matrix 
Jₙ with main diagonal a_vec = [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] 
and subdiagonal b_vec = [b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], 
where n is the number of Gauss quadrature nodes and 
[α₀, α₁, ..., αₙ₋₁] and [β₁, ..., βₙ] are the vectors of 
recursion coefficients of the three-term recurrence relation. 
The computation of these eigenvalues and eigenvectors is 
carried out using Double64 arithmetic and the package 
GenericLinearAlgebra.jl
"""
function stieltjes_step2_fn(n::Integer, μ₀, a_vec, b_vec, a, b)
  @assert n ≥ 2
  a_vec_converted = convert(Vector{Double64}, a_vec)
  b_vec_converted = convert(Vector{Double64}, b_vec)
  μ₀_converted = materialize_scalar_spec_fn(Double64, μ₀)
  Jₙ = SymTridiagonal(a_vec_converted, b_vec_converted)
  eigenvalues, eigenvectors = eigen(Jₙ)
  nodes = eigenvalues
  # println("nodes[1] = ",nodes[1] )
  # Materialize the two support endpoints separately. Rebuilding [a, b] can
  # force mixed inputs like BigFloat(0) and Inf to promote before conversion.
  T_a = materialize_scalar_spec_fn(Double64, a)
  T_b = materialize_scalar_spec_fn(Double64, b)
  @assert T_a ≤ nodes[1]
  @assert nodes[n] ≤ T_b
  weights = zeros(Double64,n)
  for j in 1:n
    weights[j] = eigenvectors[1, j]^2 / dot(eigenvectors[:,j], eigenvectors[:,j])
  end
  weights = μ₀_converted * weights
  [nodes, weights]
end


"""
nodes, weights = stieltjes_custom_gauss_quad_all_fn(n, μ₀, a_vec, b_vec, a, b, upto_n) 

This function computes the custom-made Gauss quadrature nodes and 
weights, for n nodes and the weight function specified by the inputs 
μ₀, a_vec, b_vec, a and b. This function is used only when n ≥ 2.      
The Boolean variable upto_n takes the value true if the Gauss quadrature 
rules are computed for number of nodes q from 1 up and including n.
By default, this variable is false, so that only the Gauss quadrature 
rule for number of nodes n is computed.
"""
function stieltjes_custom_gauss_quad_all_fn(n::Integer, μ₀, a_vec, b_vec, 
    a, b, upto_n::Bool=false)
  @assert n ≥ 2
  nodes, weights = stieltjes_step2_fn(n, μ₀, a_vec, b_vec, a, b)
  if upto_n == false
    # The precision of BigFloat is set globally.
    # Consequently, this precision is reset to
    # its default value before exiting this
    # function. 
    setprecision(BigFloat, 256, base=2)
    return([nodes, weights])
  end
  nodes_upto_n = Any[]
  weights_upto_n = Any[]
  #  a_vec = [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] 
  nodes_1 = convert(Double64, a_vec[1])      
  weights_1 = materialize_scalar_spec_fn(Double64, μ₀)          
  push!(nodes_upto_n, nodes_1)
  push!(weights_upto_n, weights_1)
  for q in 2:(n-1)
    a_vec_q = a_vec[1:q]
    b_vec_q = b_vec[1:q]
    nodes_q, weights_q =  stieltjes_step2_fn(q, μ₀, a_vec_q, b_vec_q, a, b)
    push!(nodes_upto_n, nodes_q)
    push!(weights_upto_n, weights_q)
  end
  push!(nodes_upto_n, nodes)
  push!(weights_upto_n, weights)
  # The precision of BigFloat is set globally.
  # Consequently, this precision is reset to
  # its default value before exiting this
  # function. 
  setprecision(BigFloat, 256, base=2)

  [nodes_upto_n, weights_upto_n]
end
