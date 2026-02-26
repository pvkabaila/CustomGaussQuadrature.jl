# The latest version of my Julia package CustomGaussQuadrature is in the folder:
# \RESEARCH - NUMERICAL METHODS\QUADRATURE\Custom GAUSS\CustomGaussQuadrature - Julia package\
# CustomGaussQuadrature\

# This script tests the version on my computer, not the
# version in the Julia General Registry.

# To run the code in this package, the first step in VS code is 
# File > Open Folder... > open the above folder

# Use CTRL + Shift + P to get to the Command Palette in VS code and choose
# Julia: Start REPL
# to start the Julia REPL.

# Use 
# Go to package mode by typing ]
# Activate the package environment by running the 
# activate . command
# to activate the package. Exit package mode using Ctrl + C

# Then use the following command to run this script:
# julia> include("test/test_stjieltjes_lnf_stored_scr.jl")

# Copy the output at the REPL into a 
# text document. Do NOT execute this 
# code line-by-line in the REPL, as this introduces
# ugly extra spaces.

# Have a look at the README.md file by opening this file in VS code and then using
# CTRL + K      V 
# to get a preview of this file.

println("Test of the version of the CustomGaussQuadrature package")
println("on my computer, not the version in the Julia General Registry.")

using CustomGaussQuadrature
# This script tests the version on my computer, not the
# version in the Julia General Registry. See path.
pathof(CustomGaussQuadrature)

using Printf
using Plots
using GaussQuadrature
using SpecialFunctions

println("\n", "------------------------------------------------------")
#-------------------------------
# scaled chi pdf weight function
#-------------------------------

m = 160;
which_f = ["scaled chi pdf", [0,Inf], m];
n = 33;

println("which_f = ", which_f)
println("n = ", n, "\n")

#  @__DIR__ gives the folder of the file currently being executed
#  .. = parent folder (package root),
#  then into src,
#  then to stjieltjes_lnf_stored_scr.jl

include(joinpath(@__DIR__, "..", "src", "stjieltjes_lnf_stored_scr.jl"));

#-----------------------------------------------------------------
# Compute the custom Gauss quadrature nodes and weights
# using the moment determinant method
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

# Compare the nodes and weights computed by the Stjieltjes
# procedure with the nodes and weights computed using the
# moment determinant method
diff_nodes = stjieltjes_nodes - nodes;
rel_diff_weights = (stjieltjes_weights - weights) ./ weights;

diff_nodes = convert(Vector{Float64}, diff_nodes);
rel_diff_weights = convert(Vector{Float64}, rel_diff_weights);

println("      stjieltjes_nodes - nodes    (stjieltjes_weights - weights)./weights")
for i in 1:lastindex(diff_nodes)
	@printf "%2d     " i
	@printf "%.16e     " diff_nodes[i]
    @printf "%.16e  \n" rel_diff_weights[i]
end

println("\n", "------------------------------------------------------")
#-------------------------------
# Hermite weight function
#-------------------------------

which_f = ["Hermite", [-Inf,Inf]];
n = 10;

println("which_f = ", which_f)
println("n = ", n, "\n")

#  @__DIR__   gives the folder of the file currently being executed
#  ..         takes us up one directory to the parent directory of this folder
#  src        takes us into the src folder
#  then to stjieltjes_lnf_stored_scr.jl

#  The following command results in values for 
#  stjieltjes_nodes and stjieltjes_weights

include(joinpath(@__DIR__, "..", "src", "stjieltjes_lnf_stored_scr.jl"));

#-----------------------------------------------------------------
# Compute the custom Gauss quadrature nodes and weights
# using the moment determinant method
moment_fn = moment_weibull_pdf_fn
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

# Compare the nodes and weights computed by the Stjieltjes
# procedure with the nodes and weights computed using the
# moment determinant method
diff_nodes = stjieltjes_nodes - nodes;
rel_diff_weights = (stjieltjes_weights - weights) ./ weights;

diff_nodes = convert(Vector{Float64}, diff_nodes);
rel_diff_weights = convert(Vector{Float64}, rel_diff_weights);

println("      stjieltjes_nodes - nodes    (stjieltjes_weights - weights)./weights")
for i in 1:lastindex(diff_nodes)
	@printf "%2d     " i
	@printf "%.16e     " diff_nodes[i]
    @printf "%.16e  \n" rel_diff_weights[i]
end



