# This Julia script is used for computing the 
# custom Gauss quadrature nodes and weights
# using the Stjieltes procedure.
# This Julia script is NOT part of the module 
# CustomGaussQuadrature.jl. Consequently, it
# will NOT run when the package CustomGaussQuadrature
# is loaded. This is in accordance with the following
# general advice:

# In a Julia package, if you have a script that doesn't 
# define any functions and you want to prevent it from 
# running when the package is loaded, you should place 
# it in the src directory but do not include it in the 
# main file (typically src/MyPackage.jl).

# Recommended Approach
# Create a Separate File: Place your script in a separate 
# file within the src directory. For example, you could 
# name it src/my_script.jl.
# Do Not Include in Main File: Ensure that you do not 
# include or use the include("my_script.jl") command in 
# your package's main file (src/MyPackage.jl). This way, 
# the script won't run when the package is loaded.
# Example Structure

#  MyPackage/
#  │
#  ├── src/
#  │   ├── MyPackage.jl       # Main package file
#  │   └── my_script.jl       # Your script
#  │
#  ├── test/
#  │   └── runtests.jl        # Tests for the package
#  │
#  └── Project.toml           # Package configuration
#  Summary
#  By structuring your package this way, the script will only
#  run when explicitly called from other scripts or the Julia REPL, 
#  ensuring it does not execute upon package loading.

if which_f[1] == "scaled chi pdf"
    moment_fn = moment_stored_fn;
    m = which_f[3];
    @assert m > 0
    lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
    T = BigFloat;
    μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
    l_endpt, u_endpt = which_f[2];
    T_l_endpt = parse(T, string(l_endpt));
    T_u_endpt = parse(T, string(u_endpt)); 
elseif which_f[1] == "Hermite"
    moment_fn = moment_stored_fn;
    lnf_fn = x -> lnf_hermite_fn(x);
    T = BigFloat;
    μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
    l_endpt, u_endpt = which_f[2];
    T_l_endpt = parse(T, string(l_endpt));
    T_u_endpt = parse(T, string(u_endpt)); 
elseif which_f[1] == "Generalized Laguerre"
    moment_fn = moment_stored_fn;
    α_GGL = which_f[3]
    @assert α_GGL > -1
     T = BigFloat;
    T_α_GGL = parse(T, string(α_GGL))
    lnf_fn = x -> lnf_laguerre_fn(x, T_α_GGL);
    μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
    l_endpt, u_endpt = which_f[2];
    T_l_endpt = parse(T, string(l_endpt));
    T_u_endpt = parse(T, string(u_endpt)); 
elseif which_f[1] ==  "chemistry example"
    moment_fn = moment_stored_fn;
    lnf_fn = x -> lnf_chemistry_fn(T, x);
    T = BigFloat;
    μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
    l_endpt, u_endpt = which_f[2];
    T_l_endpt = parse(T, string(l_endpt));
    T_u_endpt = parse(T, string(u_endpt)); 
else
    # DomainError means that the argument to a function
    # or constructor does not lie in the valid domain
    throw(DomainError(which_f, "invalid argument"))
end

a = T_l_endpt;
b = T_u_endpt;
stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);
