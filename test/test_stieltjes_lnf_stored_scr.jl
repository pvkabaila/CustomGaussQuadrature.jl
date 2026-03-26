# This script tests the high-level script stieltjes_lnf_stored_scr.jl
# which resides in 
# EITHER 
# (1) The local development version of my Julia package 
#     CustomGaussQuadrature, which is in the folder:
#     \RESEARCH - NUMERICAL METHODS\QUADRATURE\Custom GAUSS\
#     CustomGaussQuadrature - Julia package\CustomGaussQuadrature\
# OR
# (2) The version of my Julia package CustomGaussQuadrature in the 
#     Julia General Registry.
#
# Choose path (1) or (2) below and follow the corresponding
# instructions.
#
#------------------------------------------------------------
# How to run this script — Path (1): Local development version
#------------------------------------------------------------
#
# Step 1 — Open the folder for the local development version in VS Code:
#     File > Open Folder... > open the folder for the local 
#     development version of my Julia package CustomGaussQuadrature, 
#     which is in the folder:
#     \RESEARCH - NUMERICAL METHODS\QUADRATURE\Custom GAUSS\
#     CustomGaussQuadrature - Julia package\CustomGaussQuadrature\
#
# Step 2 — Start the Julia REPL:
#     CTRL + Shift + P  >  Julia: Start REPL
#
# Step 3 — Activate the local package:
#     julia> ]
#     (@v1.xx) pkg> activate .
#     (CustomGaussQuadrature) pkg>    (press Ctrl + C to exit pkg mode)
#
# Step 4 — Run the script:
#     julia> using CustomGaussQuadrature
#     julia> include("test/test_stieltjes_lnf_stored_scr.jl")
#
#------------------------------------------------------------
# How to run this script — Path (2a): Julia General Registry
# version, when the folder for the local development version 
# of the CustomGaussQuadrature package is already open in VS Code, 
# but is NOT activated
#------------------------------------------------------------
#
# Step 1 — Start the Julia REPL:
#     CTRL + Shift + P  >  Julia: Start REPL
#
# Step 2 — Activate a temporary environment, install the package,
#           and load it:
#     julia> using Pkg
#     julia> Pkg.activate(mktempdir())
#     julia> Pkg.add("CustomGaussQuadrature")
#     julia> Pkg.add("Plots")
#     julia> Pkg.add("GaussQuadrature")
#     julia> Pkg.add("SpecialFunctions")
#     julia> using CustomGaussQuadrature
#
#   A temporary environment is necessary because Julia's package
#   manager does not allow Pkg.add("CustomGaussQuadrature") when
#   the active project already has the same name or UUID.
#
# Step 4 — Run the script:
#     julia> pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
#     julia> include(joinpath(pkg_dir, "test", "test_stieltjes_lnf_stored_scr.jl"))
#
#------------------------------------------------------------
# How to run this script — Path (2b): Julia General Registry
# version, when the folder for the local development version 
# of the CustomGaussQuadrature package is NOT open in VS Code
#------------------------------------------------------------
#
# Step 1 — Start the Julia REPL:
#     CTRL + Shift + P  >  Julia: Start REPL
#
# Step 2 — Install the package and load it:
#     julia> using Pkg
#     julia> Pkg.add("CustomGaussQuadrature")
#     julia> Pkg.add("GaussQuadrature")
#     julia> Pkg.add("SpecialFunctions")
#     julia> using CustomGaussQuadrature
#
# Step 4 — Run the script:
#     julia> pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
#     julia> include(joinpath(pkg_dir, "test", "test_stieltjes_lnf_stored_scr.jl"))
#
#------------------------------------------------------------
# DETAILS
#
# pathof(CustomGaussQuadrature) returns the absolute path to the
# main source file CustomGaussQuadrature.jl.  Its value depends on
# which version was loaded:
#
# (a) Local version (] activate .) — returns e.g.
#     C:\Users\pkaba\...\CustomGaussQuadrature\src\CustomGaussQuadrature.jl
#     This is the source code you are editing on your computer.
#
# (b) Registry version (Pkg.add) — returns e.g.
#     C:\Users\pkaba\.julia\packages\CustomGaussQuadrature\Ab1Cd\src\CustomGaussQuadrature.jl
#     The "Ab1Cd" part is a short hash that Julia uses to identify
#     the installed version.  The source code lives inside .julia
#     and is managed by the package manager.
#
# In both cases, dirname(dirname(pathof(...))) gives the package
# root, and the include calls in this script resolve to the correct
# files.
#------------------------------------------------------------

# Copy the output at the REPL into a 
# text document. Do NOT execute this 
# code line-by-line in the REPL, as this introduces
# ugly extra spaces.

# Have a look at the README.md file by opening this file in VS code and then using
# CTRL + K      V 
# to get a preview of this file.

# The following script is used for plotting cdf's and empirical
# distribution functions as well as other functions useful for
# code development and journal article writing.
include("utilities_scr.jl")

# Both the inbuilt package Dates and the package Printf have a format command
date_time_now = now()
println("\n", Dates.Date(date_time_now), 
"   ", Dates.format(date_time_now, "HH:MM"))    

# pkg_dir is the package root directory, regardless of whether the
# package was loaded from a local path or from the registry.
pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
println("\n", "The package root is:")
println_wrap(pkg_dir)
print("\n")

println("The simplest way to compare the output from this test,")
println("stored as a .txt ('plain text') file, with a previous output")
println("from this test is as follows. Open both .txt files")
println("in VS Code. Then on the EXPLORER menu on the left")
println("right click on the previous output > Select for Compare.")
println("Then right click on the .txt file for the latest output >")
println("Compare with Selected. Both .txt files are compared.")
println("This does not compare digits in a number numerically,")
println("by first rounding to the appropriate number of digits")
println("and then comparing, it simply compares the difference")
println("between characters. Nonetheless, VS code diff has the")
println("virtue of being automated.")
println("Unfortunately, VS code diff will pick up differences that")
println("are invisible to me, for example trailing whitespaces.")
println("These trailing whitespaces can be rendered by opening")
println("the Command Palette with Ctrl+Shift+P, typing") 
println("'Toggle Render Whitespace' and pressing Enter.")
println("Also, the fact that a difference is due to a trailing")
println("whitespace can actually be perceived by careful")
println("examination of the diff view (check the dotted lines).")
println("However, the best approach is to first use the following")
println("Powershell script on both .txt files:")
println(".\\remove_trailing_whitespace.ps1 -FilePath 'filename.txt'")
println("written by Claude Sonnet 4.6.", "\n")

println("2026 3 25 I used Claude Opus 4.6 to check the safety of installing")
println("the EXTENSION Persistent Text Marker and it reported that no") 
println("dangerous, suspicious or questionable code was found.") 
println("Unfortunately, this highlighting does not work the way")
println("that I want it to for highlighting part of a number")
println("since this highlights ALL occurences in the file of")
println("that part of the number. This cannot be changed.")
println("The only non-confusing way that this highlighter")
println("can be used is if I add comments which I highlight")
println("in Yellow. I do this by right clicking on the chosen")
println("and then choosing the colour Yellow.", "\n")

println("To use VS code diff to aid in highlighting differences")
println("in the original .txt file, I open up this file in a normal")
println("editor split alongside the diff view, so I can see  diff")
println("view, while applying the highlights in the other pane.")
println("To go to the same line in both diff view & regular editor")
println("use Ctrl+G, type the line number and press Enter - works")
println("in both the diff view and the regular editor tab to jump") 
println("quickly to the same position.", "\n")

using Printf
using GaussQuadrature
using SpecialFunctions
using CustomGaussQuadrature: μ_offsetvec_fn

println("\n", "------------------------------------------------------")
#-------------------------------
# scaled chi pdf weight function
#-------------------------------

m = 160;
which_f = ["scaled chi pdf", [0,Inf], m];
n = 33;

println("which_f = ", which_f)
println("n = ", n, "\n")

include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

println("r =", r, "\n")

#-----------------------------------------------------------------
# Compute the custom Gauss quadrature nodes and weights
# using the moment determinant method
moment_fn = moment_stored_fn
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

# Compare the nodes and weights computed by the Stieltjes
# procedure with the nodes and weights computed using the
# moment determinant method
println("Compare the nodes and weights computed by the Stieltjes")
println("procedure with the nodes and weights computed using the")
println("moment determinant method.", "\n")
diff_nodes = stieltjes_nodes - nodes;
rel_diff_weights = (stieltjes_weights - weights) ./ weights;

diff_nodes = convert(Vector{Float64}, diff_nodes);
rel_diff_weights = convert(Vector{Float64}, rel_diff_weights);

println("      stieltjes_nodes - nodes    (stieltjes_weights - weights)./weights")
for i in 1:lastindex(diff_nodes)
	@printf "%2d     " i
	@printf "%.16e     " diff_nodes[i]
    @printf "%.16e  \n" rel_diff_weights[i]
end

print("\n")
println("maximum(abs.(stieltjes_nodes - nodes)) = ", 
convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((stieltjes_weights - weights) ./ weights)) = ", 
convert(Float64, maximum(abs.(rel_diff_weights))))
print("\n")

println("\n", "------------------------------------------------------")
#-------------------------------
# Hermite weight function
#-------------------------------
# Compute the custom Gauss quadrature nodes and weights
# using the command hermite from the GaussQuadrature package
# with arithmetic type BigFloat 

which_f = ["Hermite", [-Inf,Inf]];
n = 10;

println("which_f = ", which_f)
println("n = ", n, "\n")

include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

println("r =", r, "\n")

println("setprecision(BigFloat, 256, base=2) ", "\n")
setprecision(BigFloat, 256, base=2);

println("Use the command hermite from the GaussQuadrature package:")
println("nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);")
nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
print("\n")

#-----------------------------------------------------------------
# Compare the nodes and weights computed by the Stieltjes
# procedure with the nodes and weights computed using the
# command hermite from the GaussQuadrature package 
# with arithmetic type BigFloat
diff_nodes = stieltjes_nodes - nodes_BigFloat;
rel_diff_weights = (stieltjes_weights - weights_BigFloat) ./ weights_BigFloat;

diff_nodes = convert(Vector{Float64}, diff_nodes);
rel_diff_weights = convert(Vector{Float64}, rel_diff_weights);

println("Comparison with nodes and weights that are, to Double64 precision, exact.", "\n")
println("      stieltjes_nodes - nodes    (stieltjes_weights - weights)./weights")
for i in 1:lastindex(diff_nodes)
	@printf "%2d     " i
	@printf "%.16e     " diff_nodes[i]
    @printf "%.16e  \n" rel_diff_weights[i]
end

print("\n")
println("maximum(abs.(stieltjes_nodes - nodes)) = ", 
convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((stieltjes_weights - weights) ./ weights)) = ", 
convert(Float64, maximum(abs.(rel_diff_weights))))

