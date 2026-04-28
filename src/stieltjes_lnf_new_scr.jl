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

a, b = which_f[2];
# Keep the support endpoints in their input representation.
# The Stieltjes driver chooses the working type T internally, and the
# support-quadrature helpers convert a and b to that local T when needed.
# This prevents the stored type of the endpoints from forcing BigFloat or
# Double64 arithmetic before the driver has chosen which branch to test.

μ₀ = mu0;

# This script keeps lnf_typed_fn and which_f separate and passes both to the
# driver. The driver then chooses T and creates the final one-argument lnf_fn.
# In closure language, the driver calls a function factory that takes T and
# returns x -> lnf_typed_fn(T, which_f, x). That returned function has x as
# its only explicit argument, but it still has access to T and which_f from
# the scope in which the closure was created.
if @isdefined k_max
    stieltjes_final_user_fn = (node_count) ->
    CustomGaussQuadrature.stieltjes_a_vec_b_vec_final_fn(node_count, μ₀, lnf_typed_fn, which_f, a, b, 7, k_max)
else
    stieltjes_final_user_fn = (node_count) ->
    CustomGaussQuadrature.stieltjes_a_vec_b_vec_final_fn(node_count, μ₀, lnf_typed_fn, which_f, a, b)
end

@assert n ≥ 1
if n == 1
    stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, stieltjes_r = 
    stieltjes_final_user_fn(2);
    stieltjes_nodes = convert(Double64, stieltjes_a_vec[1])
    stieltjes_weights = convert(Double64, 1)
else
    stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, stieltjes_r = 
    stieltjes_final_user_fn(n);

    stieltjes_nodes, stieltjes_weights = 
    stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b);
end