# This Julia script is part of the module 
# CustomGaussQuadrature.jl
#
# This script consists of some utility functions for the
# Julia package CustomGaussQuadrature.jl

"""
plot_cdf_discrete_rv_fn(x_vec, prob_vec, x_lo, x_hi)

Returns a graph of the cumulative distribution function (cdf) of a 
discrete random variable with point probability masses prob_vec
at locations x_vec. The horizontal axis of this graph is for
x between x_lo and x_hi.
"""
function plot_cdf_discrete_rv_fn(x_vec, prob_vec, x_lo, x_hi)
    cdf_at_point_masses = cumsum(prob_vec)
    plot(x_vec, cdf_at_point_masses, seriestype=:scatter, mc=:blue, ms=1,
    legend=false)
    xlims!(x_lo, x_hi)
    ylims!(0.0, 1.05)
    xlabel!("x")
    ylabel!("cdf(x)")
    hline!([1.0], lc=:black)
    for i = 1:(length(x_vec)-2)
        x_i = x_vec[i]
        x_ip1 = x_vec[i+1]
        cdf_val = cdf_at_point_masses[i]
        plot!([x_i, x_ip1], [cdf_val, cdf_val], legend=false, lc=:black)
    end
    cdf_val = cdf_at_point_masses[length(x_vec)]
    plot!([x_vec[length(x_vec)-1], x_vec[length(x_vec)]], [cdf_val, cdf_val], 
    legend=false, lc=:black)
    
end


"""
scaled_chi_cdf = scaled_chi_cdf_fn(x, m)

Returns the cumulative distribution function (cdf)
of a random variable with the the scaled chi pdf
with m degrees of freedom, evaluated at x, using 
the known cdf of a χ² distribution with m degrees
of freedom.
"""
function scaled_chi_cdf_fn(x, m::Integer)
    cdf(Chisq(m), m * x^2)
end


"""
weibull_cdf = weibull_cdf_fn(x, k)

Returns the cumulative distribution function (cdf)
of a random variable with the the a Weibull distribution
with scale parameter λ=1 and shape parameter k (k > 0).
"""
function weibull_cdf_fn(x, k)
    @assert k > 0
    if x < 0 
        return(0)
    end
    1 - exp(-x^k)
end



"""
integrand_transf_integral = integrand_transf_integral_fn(T, y, k, f_fn, a, b)

Returns ϕᵏ(y) f(ϕ(y)) ϕ′(y), where f is a specified nonnegative 
weight function with support interval [0, ∞). The function returned 
is the integrand of the right-hand side of
    ∞               1
    ∫ xᵏ f(x) dx  = ∫ ϕᵏ(y) f(ϕ(y)) ϕ′(y) dy, 
    0              -1
where k is a specified positive integer. 
"""
function integrand_transf_integral_fn(T, y, k::Integer, f_fn::Function, a, b)
    @assert convert(T, -1) ≤ y < convert(T, 1)
    ϕ_val = ϕ_fn(T, y, a, b) 
    ϕ_val^k * f_fn(ϕ_val) * ϕ′_fn(T, y, a, b) 
end







"""
nodes_upto_n_Float64, weights_upto_n_Float64 =
nodesweights_upton_to_Float64_fn(nodes_upto_n, weights_upto_n, n)
"""
function nodesweights_upton_to_Float64_fn(nodes_upto_n, weights_upto_n, n)
    nodes_upto_n_Float64 = Vector{Vector}(undef, n)
    nodes_upto_n_Float64[1] = [convert(Float64, nodes_upto_n[1])]
    for q in 2:n
        nodes_upto_n_Float64[q] = convert(Vector{Float64}, nodes_upto_n[q])
    end
    weights_upto_n_Float64 = Vector{Vector}(undef, n)
    weights_upto_n_Float64[1] = [convert(Float64, weights_upto_n[1])]
    for q in 2:n
        weights_upto_n_Float64[q] = convert(Vector{Float64}, weights_upto_n[q])
    end
    
    [nodes_upto_n_Float64, weights_upto_n_Float64]
end


