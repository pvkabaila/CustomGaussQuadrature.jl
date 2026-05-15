# choice_of_r_scr.jl
# This script is for investigation, not for automated testing.

using Pkg
# Activate the package root, i.e. the directory one level above this test
# script that contains Project.toml for CustomGaussQuadrature.
#
# After this activation, `using CustomGaussQuadrature` is resolved in the
# active project environment defined by that Project.toml, not by whatever
# package version may be installed in the default global environment.
# Because the active project is the CustomGaussQuadrature package itself,
# Julia loads the module from this checkout's src/CustomGaussQuadrature.jl.
# In other words, the package name is satisfied by the currently activated
# local project, so the registered/installed copy is not the one that gets
# imported for this session.
Pkg.activate(joinpath(@__DIR__, ".."))
using CustomGaussQuadrature
# μ_offsetvec_fn and stieltjes_a_vec_b_vec_final_fn are internal functions 
# of CustomGaussQuadrature; the qualified import below makes these explicitly 
# available in this script.
using CustomGaussQuadrature: μ_offsetvec_fn, stieltjes_a_vec_b_vec_final_fn
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
which_f = ["scaled chi pdf", [0,Inf], m];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3]);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
    r_vec[n - 1] = r
    println("number of Gauss quadrature nodes n = ", n, ",   r = ", r, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("scaled chi pdf weight function, with m = $m, 
β_hat = $β_hat", titlefontsize=10)
savefig("scaled_chi_pdf_m=$m.pdf")





println("chemistry example weight function")
which_f = ["chemistry example", [0, Inf]];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_chemistry_fn(T, x);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
    r_vec[n - 1] = r
    println("number of Gauss quadrature nodes n = ", n, ",   r = ", r, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("chemistry weight function, 
β_hat = $β_hat", titlefontsize=10)
savefig("chemistry.pdf")



which_f = ["hermite", [-Inf, Inf]];
a = -Inf;
b = Inf;
lnf_typed_fn = (T, which_f, x) -> lnf_hermite_fn(T, x);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
    r_vec[n - 1] = r
    println("number of Gauss quadrature nodes n = ", n, ",   r = ", r, "\n")
end

scatter(n_vec, r_vec, legend=false)
xlabel!("n")
ylabel!("r")

β_hat = add_ls_stline_fn(n_vec, r_vec);
title!("Hermite weight function, 
β_hat = $β_hat", titlefontsize=10)
savefig("hermite.pdf")

#-------------------------------------------------
# Examine the effect of using an explicit offset and a larger j_max
# while keeping the driver-owned T/closure path.

m = 1;
which_f = ["scaled chi pdf", [0,Inf], m];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3]);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
j_max = 400;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=offset, j_max=j_max);
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

println("chemistry example weight function")
which_f = ["chemistry example", [0, Inf]];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_chemistry_fn(T, x);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
j_max = 400;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=offset, j_max=j_max);
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



which_f = ["hermite", [-Inf, Inf]];
a = -Inf;
b = Inf;
lnf_typed_fn = (T, which_f, x) -> lnf_hermite_fn(T, x);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
j_max = 400;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=offset, j_max=j_max);
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


m = 5;
which_f = ["scaled chi pdf", [0,Inf], m];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3]);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
j_max = 400;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=offset, j_max=j_max);
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


m = 160;
which_f = ["scaled chi pdf", [0,Inf], m];
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> ln_scaled_chi_pdf_fn(T, x, which_f[3]);
T = BigFloat;
moment_fn = moment_stored_fn;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

n_upper = 40;
r_vec = zeros(n_upper - 1);
n_vec = 2:n_upper;
offset = 7;
println("offset = ", offset)
j_max = 400;
for n in n_vec
    @time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
    stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b; offset=offset, j_max=j_max);
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


