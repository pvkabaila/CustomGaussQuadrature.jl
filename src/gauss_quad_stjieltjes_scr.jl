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
# Stjieltjes procedure, where the inner product 
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
a_vec, b_vec = stjieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_support, lnweights_support)

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
function stjieltjes_a_vec_b_vec_nonan_fn(T::Type, n::Integer, μ₀, nodes_support, lnweights_support)
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



"""
a_vec, b_vec, nbits = stjieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_fn, a, b, r)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and 
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], where n is the number 
of Gauss quadrature nodes (n ≥ 2). αₖ and βₖ are defined by 
  αₖ = (x πₖ, πₖ)ₐ / (πₖ, πₖ)ₐ and  βₖ = (πₖ, πₖ)ₐ / (πₖ₋₁, πₖ₋₁)ₐ 
where (u,v)ₐ denotes the discrete approximation               
  r                               b      
  ∑ wᵢ u(xᵢ) v(xᵢ)   to   (u,v) = ∫ u(x) v(x) f(x) dx.           
 i=1                              a                        
The algorithm chooses the precision nbits of BigFloat 
arithmetic that is used in all of the computations, so 
that there is some assurance that the vectors [a₁, ..., aₙ] 
and [b₁, ..., bₙ₋₁] are computed to sufficient accuracy,
for the specified value of r.
"""
function stjieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_fn, which_f, r)
 
    setprecision(BigFloat, 256, base=2)
    T = BigFloat
    nbits = 256
    # println("T = ", T, ",  nbits = ", precision(BigFloat))
    # @time "Compute nodes_256 & lnweights_256" 
    nodes_256, lnweights_256 = 
    nodes_lnweights_support_fn(T, lnf_fn::Function, which_f, r::Integer)
    # @time "Compute a_vec_256 & b_vec_256" 
    a_vec_256, b_vec_256 = 
    stjieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_256, lnweights_256)

    epsilon = 1.0e-22

    T = Double64
    # println("T = ", T)
    # @time "Compute nodes_106 & lnweights_support_106" 
    nodes_106, lnweights_106 = 
    nodes_lnweights_support_fn(T, lnf_fn::Function, which_f, r::Integer)
    # @time "Compute a_vec & b_vec" 
    a_vec_106, b_vec_106 = 
    stjieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_106, lnweights_106)
    if any(isnan, a_vec_106)==false && any(isnan, b_vec_106)==false
        abs_error_a = convert(Vector{Float64}, abs.(a_vec_106 - a_vec_256))
        max_abs_error_a = maximum(abs_error_a)
        # println("max_abs_error_a = ", max_abs_error_a)
        rel_error_b = convert(Vector{Float64}, 
        (abs.((b_vec_106 - b_vec_256) ./ b_vec_256)))
        max_rel_error_b = maximum(rel_error_b)
        # println("max_rel_error_b = ", max_rel_error_b)
        if max_abs_error_a < epsilon && max_rel_error_b < epsilon
            return([a_vec_256, b_vec_256, nbits])
        end
    end

    setprecision(BigFloat, 224, base=2)
    T = BigFloat;
    # println("T = ", T, ",  nbits = ", precision(BigFloat))
    # @time "Compute nodes_224 & lnweights_224" 
    nodes_224, lnweights_224 = 
    nodes_lnweights_support_fn(T, lnf_fn::Function, which_f, r::Integer)
    # @time "Compute a_vec & b_vec" 
    a_vec_224, b_vec_224 = 
    stjieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_224, lnweights_224) 
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_224 - a_vec_256))
    max_abs_error_a = maximum(abs_error_a)
    # println("max_abs_error_a = ", max_abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_224 - b_vec_256) ./ b_vec_256)))
    max_rel_error_b = maximum(rel_error_b)
    # println("max_rel_error_b = ", max_rel_error_b)
    if max_abs_error_a < epsilon && max_rel_error_b < epsilon
        return([a_vec_256, b_vec_256, nbits])
    end
            
    # println(" ")
    setprecision(BigFloat, 512, base=2)
    T = BigFloat
    nbits = 512
    # println("T = ", T, ",  nbits = ", precision(BigFloat))
    # @time "Compute nodes_512 & lnweights_512" 
    nodes_512 , lnweights_512  = 
    nodes_lnweights_support_fn(T, lnf_fn::Function, which_f, r::Integer)
    # @time "Compute a_vec & b_vec" 
    a_vec_512 , b_vec_512  = 
    stjieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_512 , lnweights_512 ) 
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_512 - a_vec_256))
    max_abs_error_a = maximum(abs_error_a)
    # println("max_abs_error_a = ", max_abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_512 - b_vec_256) ./ b_vec_512)))
    max_rel_error_b = maximum(rel_error_b)
    # println("max_rel_error_b = ", max_rel_error_b)
    if max_abs_error_a < epsilon && max_rel_error_b < epsilon
        return([a_vec_512, b_vec_512, nbits])
    end
    
    setprecision(BigFloat, 480, base=2)
    T = BigFloat;
    # println("T = ", T, ",  nbits = ", precision(BigFloat))
    # @time "Compute nodes_480 & lnweights_480" 
    nodes_480, lnweights_480 = 
    nodes_lnweights_support_fn(T, lnf_fn::Function, which_f, r::Integer)
    # @time "Compute a_vec & b_vec" 
    a_vec_480, b_vec_480 = 
    stjieltjes_a_vec_b_vec_nonan_fn(T, n, μ₀, nodes_480, lnweights_480) 
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_480 - a_vec_480))
    max_abs_error_a = maximum(abs_error_a)
    # println("max_abs_error_a = ", max_abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_480 - b_vec_512) ./ b_vec_512)))
    max_rel_error_b = maximum(rel_error_b)
    # println("max_rel_error_b = ", max_rel_error_b)
    if max_abs_error_a > epsilon || max_rel_error_b > epsilon
        throw(DomainError(n, " too large"))  
    end    
    setprecision(BigFloat, 256, base=2)

    [a_vec_512, b_vec_512, nbits]  
end 




"""
a_vec, b_vec, nbits, r = stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, which_f, offset=7, k_max=40)

Returns the vectors [a₁, ..., aₙ] = [α₀, ... , αₙ₋₁] and 
[b₁, ..., bₙ₋₁] = [√β₁, ... , √βₙ₋₁], where n is the number 
of Gauss quadrature nodes (n ≥ 2). αₖ and βₖ are defined by 
  αₖ = (x πₖ, πₖ) / (πₖ, πₖ) and  βₖ = (πₖ, πₖ) / (πₖ₋₁, πₖ₋₁) 
where              
            b      
    (u,v) = ∫ u(x) v(x) f(x) dx.           
            a                        
The algorithm chooses the precision nbits of BigFloat 
arithmetic that is used in all of the computations and the 
value of r that specifies the discrete approximation               
  r                               b      
  ∑ wᵢ u(xᵢ) v(xᵢ)   to   (u,v) = ∫ u(x) v(x) f(x) dx,           
 i=1                              a                        
so that there is some assurance that the vectors [a₁, ..., aₙ] 
and [b₁, ..., bₙ₋₁] are computed to sufficient accuracy for use
in Step 2. To choose the final value of r, which is outputted,
we consider k (offset + n) for k=3, k-4, up to a maximum 
value k_max. The default values of offset and k_max are 7
and 40, respectively. 
"""
function stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, which_f, offset=7, k_max=40)
    k = 3
    r = k * (offset + n)
    # println("r = k * (offset + n), where k = ", k)
    # @time 
    a_vec, b_vec, nbits = 
    stjieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_fn, which_f, r)
    
    k = k + 1
    r = k * (offset + n)
    # println("r = k * (offset + n), where k = ", k)
    # @time 
    a_vec_new, b_vec_new, nbits_new = 
    stjieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_fn, which_f, r)
    # println("nbits_new = ", nbits_new)
    
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_new - a_vec))
    max_abs_error_a = maximum(abs_error_a)
    # println("maximum(abs.(a_vec_new - a_vec)) = ", max_abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_new - b_vec) ./ b_vec_new)))
    max_rel_error_b = maximum(rel_error_b)
    # println("maximum(abs.(b_vec_new - b_vec) ./ b_vec_new)) = ", max_rel_error_b)
    
    epsilon=1.0e-18
    
    while (max_abs_error_a > epsilon || max_rel_error_b > epsilon) 
    a_vec, b_vec, nbits = a_vec_new, b_vec_new, nbits_new
    if k == k_max
        throw(DomainError(k, " needed value of k exceeds k_max"))
    end
    k = k + 1
    r = k * (offset + n)
    # println(" ")
    # println("r = k * (offset + n), where k = ", k)
    # @time 
    a_vec_new, b_vec_new, nbits_new = 
    stjieltjes_a_vec_b_vec_choosenbits_fn(n, μ₀, lnf_fn, which_f, r)
    # println("nbits_new = ", nbits_new)
    
    abs_error_a = convert(Vector{Float64}, abs.(a_vec_new - a_vec))
    max_abs_error_a = maximum(abs_error_a)
    # println("maximum(abs.(a_vec_new - a_vec)) = ", max_abs_error_a)
    rel_error_b = convert(Vector{Float64}, 
    (abs.((b_vec_new - b_vec) ./ b_vec_new)))
    max_rel_error_b = maximum(rel_error_b)
    # println("maximum(abs.(b_vec_new - b_vec) ./ b_vec_new)) = ", max_rel_error_b)
    
    end
    
    # The precision of BigFloat is set globally.
    # Consequently, this precision is reset to
    # its default value before exitting this
    # function. 
    setprecision(BigFloat, 256, base=2)
    # println("maximum(abs.(a_vec_new - a_vec)) = ", max_abs_error_a)
    # println("maximum(abs.(b_vec_new - b_vec) ./ b_vec_new)) = ", max_rel_error_b)
    [a_vec_new, b_vec_new, nbits_new, r]
    
end


"""
nodes, weights = stjieltjes_step2_fn(n, μ₀, a_vec, b_vec, which_f)

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
function stjieltjes_step2_fn(n::Integer, μ₀, a_vec, b_vec, which_f)
  a, b = which_f[2]
  T_a = parse(T, string(a))
  T_b = parse(T, string(b))
  @assert n ≥ 2
  a_vec_converted = convert(Vector{Double64}, a_vec)
  b_vec_converted = convert(Vector{Double64}, b_vec)
  μ₀_converted = convert(Double64, μ₀)
  Jₙ = SymTridiagonal(a_vec_converted, b_vec_converted)
  eigenvalues, eigenvectors = eigen(Jₙ)
  nodes = eigenvalues
  # println("nodes[1] = ",nodes[1] )
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
nodes, weights = stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, a_vec, b_vec, a, b, upto_n)

This function computes the custom-made Gauss quadrature nodes and 
weights, for n nodes and the weight function specified by the inputs 
μ₀, μ₁, a_vec, b_vec, a and b
The Boolean variable upto_n takes the value true if the Gauss quadrature 
rules are computed for number of nodes q from 1 up and including n.
By default, this variable is false, so that only the Gauss quadrature 
rule for number of nodes n is computed.
"""
function stjieltjes_custom_gauss_quad_all_fn(n::Integer, μ₀, μ₁, a_vec, b_vec, 
    a, b, upto_n::Bool=false)
  @assert n ≥ 1
  if n == 1
    μ₀ = convert(Double64, μ₀)
    μ₁ = convert(Double64, μ₁)  
    nodes = μ₁ / μ₀
    weights = μ₀
    return([nodes, weights])
  end
  nodes, weights = stjieltjes_step2_fn(n, μ₀, a_vec, b_vec, a, b)
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
  μ₀ = convert(Double64, μ₀)
  μ₁ = convert(Double64, μ₁) 
  nodes_1 = μ₁ / μ₀
  weights_1 = μ₀
  push!(nodes_upto_n, nodes_1)
  push!(weights_upto_n, weights_1)
  for q in 2:(n-1)
    a_vec_q = a_vec[1:q]
    b_vec_q = b_vec[1:q]
    nodes_q, weights_q =  stjieltjes_step2_fn(q, μ₀, a_vec_q, b_vec_q, a, b)
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



