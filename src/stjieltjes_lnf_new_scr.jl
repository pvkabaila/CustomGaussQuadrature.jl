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

l_endpt, u_endpt = which_f[2];
T_l_endpt = parse(T, string(l_endpt));
T_u_endpt = parse(T, string(u_endpt)); 
a = T_l_endpt;
b = T_u_endpt;

stjieltjes_a_vec, stjieltjes_b_vec, stjieltjes_nbits, r = 
stjieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);

stjieltjes_nodes, stjieltjes_weights = 
stjieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stjieltjes_a_vec, stjieltjes_b_vec, a, b);