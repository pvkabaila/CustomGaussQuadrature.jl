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


function stieltjes_a_vec_b_vec_choosenbits_core_fn(n, μ₀, make_lnf_fn::Function, a, b, r)
  setprecision(BigFloat, 256, base=2)
  T = BigFloat
  nbits = 256
  # Build lnf_fn only after choosing T, so every downstream computation
  # sees a closure whose constants and branches are created for that T.
  # The support endpoints a and b are still passed in as raw values; the
  # support-quadrature helper converts them to this local T, so endpoint
  # storage does not bias the BigFloat-versus-Double64 comparison.
  lnf_fn = make_lnf_fn(T)
  nodes_256, lnweights_256 = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_256, b_vec_256 = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_256, lnweights_256)

  epsilon = 1.0e-22

  T = Double64
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
    if max_abs_error_a < epsilon && max_rel_error_b < epsilon
      return([a_vec_256, b_vec_256, nbits])
    end
  end

  setprecision(BigFloat, 224, base=2)
  T = BigFloat;
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
  if max_abs_error_a < epsilon && max_rel_error_b < epsilon
    return([a_vec_256, b_vec_256, nbits])
  end
            
  setprecision(BigFloat, 512, base=2)
  T = BigFloat
  nbits = 512
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
  if max_abs_error_a < epsilon && max_rel_error_b < epsilon
    return([a_vec_512, b_vec_512, nbits])
  end
    
  setprecision(BigFloat, 480, base=2)
  T = BigFloat;
  lnf_fn = make_lnf_fn(T)
  nodes_480, lnweights_480 = 
  nodes_lnweights_support_fn(T, lnf_fn, a, b, r)
  a_vec_480, b_vec_480 = 
  stieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_480, lnweights_480) 
  abs_error_a = convert(Vector{Float64}, abs.(a_vec_480 - a_vec_480))
  max_abs_error_a = maximum(abs_error_a)
  rel_error_b = convert(Vector{Float64}, 
  (abs.((b_vec_480 - b_vec_512) ./ b_vec_512)))
  max_rel_error_b = maximum(rel_error_b)
  if max_abs_error_a > epsilon || max_rel_error_b > epsilon
    throw(DomainError(n, " too large"))  
  end    
  setprecision(BigFloat, 256, base=2)

  [a_vec_512, b_vec_512, nbits]
end


"""
a_vec, b_vec, nbits = stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], where n is the number
of Gauss quadrature nodes (n ≥ 2).

The input lnf_typed_fn is expected to have the form
lnf_typed_fn(T, which_f, x). The algorithm first chooses the
arithmetic type T and then constructs the one-argument closure
lnf_fn(x) = lnf_typed_fn(T, which_f, x) that is used in the
support quadrature and Stieltjes recurrence computations.
The support endpoints a and b are kept as raw values until the
support quadrature converts them to this same local T, so their
stored type does not influence which arithmetic branch is used.
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
    μ₀,
    T -> (x -> lnf_typed_fn(T, which_f, x)),
    a,
    b,
    r,
  )
end


"""
a_vec, b_vec, nbits, r = stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b, offset=7, k_max=40)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], together with the chosen
BigFloat precision nbits and the final value of r.

The input lnf_typed_fn is expected to have the form
lnf_typed_fn(T, which_f, x). For each trial value of r, the
driver chooses the arithmetic type T and then creates the
one-argument closure lnf_fn(x) = lnf_typed_fn(T, which_f, x)
before carrying out the support quadrature and Stieltjes
computations.
"""
function stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b, offset=7, k_max=40)
  @assert n ≥ 2
    k = 3
    r = k * (offset + n)
    a_vec, b_vec, nbits = 
    stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

    k = k + 1
    r = k * (offset + n)
    a_vec_new, b_vec_new, nbits_new = 
    stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

    abs_error_a = convert(Vector{Float64}, abs.(a_vec_new - a_vec))
    max_abs_error_a = maximum(abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_new - b_vec) ./ b_vec_new)))
    max_rel_error_b = maximum(rel_error_b)

    epsilon=1.0e-18

    while (max_abs_error_a > epsilon || max_rel_error_b > epsilon)
    a_vec, b_vec, nbits = a_vec_new, b_vec_new, nbits_new
    if k == k_max
      throw(DomainError(k, " needed value of k exceeds k_max"))
    end
    k = k + 1
    r = k * (offset + n)
    a_vec_new, b_vec_new, nbits_new = 
    stieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_typed_fn, which_f, a, b, r)

    abs_error_a = convert(Vector{Float64}, abs.(a_vec_new - a_vec))
    max_abs_error_a = maximum(abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_new - b_vec) ./ b_vec_new)))
    max_rel_error_b = maximum(rel_error_b)

    end

    setprecision(BigFloat, 256, base=2)
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
  μ₀_converted = convert(Double64, μ₀)
  Jₙ = SymTridiagonal(a_vec_converted, b_vec_converted)
  eigenvalues, eigenvectors = eigen(Jₙ)
  nodes = eigenvalues
  # println("nodes[1] = ",nodes[1] )
  @assert a ≤ nodes[1]
  @assert nodes[n] ≤ b
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
    # its default value before exitting this
    # function. 
    setprecision(BigFloat, 256, base=2)
    return([nodes, weights])
  end
  nodes_upto_n = Any[]
  weights_upto_n = Any[]
  #  a_vec = [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] 
  nodes_1 = convert(Double64, a_vec[1])      
  weights_1 = convert(Double64, μ₀)          
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
  # its default value before exitting this
  # function. 
  setprecision(BigFloat, 256, base=2)

  [nodes_upto_n, weights_upto_n]
end
