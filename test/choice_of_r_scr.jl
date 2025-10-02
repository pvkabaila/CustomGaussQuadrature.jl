# choice_of_r_scr.jl

using CustomGaussQuadrature
pathof(CustomGaussQuadrature)
using Plots


"""
add_ls_stline_fn(x, y)

x and y are vectors of the same length.
Add a straight line fitted by least squares
to scatter(x,y).
In Julia the QR decomposition is
built into the backslash operator \
"""
function add_ls_stline_fn(x, y)
    n = length(y)
    X = hcat(ones(n), x)
    β_hat = X\y
    y_lo = transpose([1, x[1]]) * β_hat
    y_hi = transpose([1, x[n]]) * β_hat
    plot!([x[1], x[n]], [y_lo, y_hi])
    β_hat
end

# m = 160;
# m = 5;
m = 1;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
for n in n_vec
    @time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, k = 
    stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
    r_vec[n - 1] = k*n
    println("number of Gauss quadrature nodes n = ", n, ",   k = ", k, ",   r = kn = ", k*n, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("scaled chi pdf weight function, with m = $m, 
β_hat = $β_hat", titlefontsize=10)
savefig("scaled_chi_pdf_m=$m.pdf")





lnf_fn = x -> lnf_chemistry_fn(T, x);
println("chemistry example weight function")
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
for n in n_vec
    @time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, k = 
    stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
    r_vec[n - 1] = k*n
    println("number of Gauss quadrature nodes n = ", n, ",   k = ", k, ",   r = kn = ", k*n, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("chemistry weight function, 
β_hat = $β_hat", titlefontsize=10)
savefig("chemistry.pdf")



which_f = ["Hermite", [-Inf, Inf]];
lnf_fn = lnf_hermite_fn;
a = -Inf;
b = Inf;
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
for n in n_vec
    @time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, k = 
    stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
    r_vec[n - 1] = k*n
    println("number of Gauss quadrature nodes n = ", n, ",   k = ", k, ",   r = kn = ", k*n, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("Hermite weight function, 
β_hat = $β_hat", titlefontsize=10)
savefig("hermite.pdf")

#-------------------------------------------------
# Examine the effect of changing r = kn, for k=3,...
# to r = k(offset+n), for k=3,...

# m = 160;
# m = 5;
m = 1;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
k_max = 400;
for n in n_vec
    @time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
    stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b, offset, k_max);
    r_vec[n - 1] = r
    println("number of Gauss quadrature nodes n = ", n, ",   r = ", r, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("scaled chi pdf weight function, with m = $m, offset = $offset,
β_hat = $β_hat", titlefontsize=10)
savefig("scaled_chi_pdf m=$m offset=$offset.pdf")

#----------------------------------------------

lnf_fn = x -> lnf_chemistry_fn(T, x);
println("chemistry example weight function")
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
k_max = 400;
for n in n_vec
    @time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
    stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b, offset, k_max);
    r_vec[n - 1] = r
    println("number of Gauss quadrature nodes n = ", n, ",   r = ", r, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("chemistry weight function, offset = $offset,
β_hat = $β_hat", titlefontsize=10)
savefig("chemistry offset=$offset.pdf")



which_f = ["Hermite", [-Inf, Inf]];
lnf_fn = lnf_hermite_fn;
a = -Inf;
b = Inf;
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
k_max = 400;
for n in n_vec
    @time "stjieltjes_a_vec_b_vec_final_fn" stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
    stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
    r_vec[n - 1] = r
    println("number of Gauss quadrature nodes n = ", n, ",   r = ", r, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("Hermite weight function, offset = $offset,
β_hat = $β_hat", titlefontsize=10)
savefig("hermite offset=$offset.pdf")
