# This script tests the functions for computing the Gauss
# rule using the (a) moment determinant method and (b) 
# Stieltjes procedure.
#
# Background to this file
#
# Apart from runtests.jl, this was the first substantial test script
# written for CustomGaussQuadrature.
#
# It was originally used as a development script for checking low-level
# and high-level computations for the moment determinant method.
# Later it was expanded to include checks for the Stieltjes procedure.
#
# Over time this file also came to be used as an informal manual
# regression-check script: after changing the package, I rerun this
# script and compare its output with previously saved output files.
#
# Consequently, this file served as a development script and is a
# practical record of how the package was checked during its
# development.
#
# There are 3 paths for running this script, which differ by:
# (i)  which test script is used (local or registry), and
# (ii) which version of CustomGaussQuadrature is under test
#      (local or registry).
#
#   Path (1):  Test script LOCAL    |  Package under test LOCAL
#   Path (2a): Test script LOCAL    |  Package under test REGISTRY
#   Path (2b): Test script REGISTRY |  Package under test REGISTRY
#
# Choose a path below and follow the corresponding instructions.
#
#------------------------------------------------------------
# How to run this script — Path (1): Local development version
#    => Test script: LOCAL  |  Package under test: LOCAL
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
#     Ctrl+Shift+P  (to bring up Command Palette)>  Julia: Start REPL
#
# Step 3 — Activate the local package:
#     julia> ]
#     (@v1.xx) pkg> activate .
#     (CustomGaussQuadrature) pkg>    (press Ctrl + C to exit pkg mode)
#
# Step 4 — Run the script:
#     julia> using CustomGaussQuadrature
#     julia> include("test/test_CustomGaussQuadrature_scr.jl")
#
#------------------------------------------------------------
# How to run this script — Path (2a): Julia General Registry
# version, when the folder for the local development version 
# of the CustomGaussQuadrature package is already open in VS Code, 
# but is NOT activated
#    => Test script: LOCAL  |  Package under test: REGISTRY
#------------------------------------------------------------
#
# Step 1 — Start the Julia REPL:
#     Ctrl+Shift+P  (to bring up Command Palette)>  Julia: Start REPL
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
#   If I want to, I can check the package version using:
#   julia> pkgversion(CustomGaussQuadrature)
#
# Step 3 — Run the script:
#     julia> include("test/test_CustomGaussQuadrature_scr.jl")
#
#   This works because the VS Code Julia REPL's working directory
#   is the local project root (the folder open in VS Code).
#   The local test script is used, but the registry version of
#   the package is under test (loaded in Step 2).
#
#------------------------------------------------------------
# How to run this script — Path (2b): Julia General Registry
# version, when the folder for the local development version 
# of the CustomGaussQuadrature package is NOT open in VS Code
#    => Test script: REGISTRY  |  Package under test: REGISTRY
#------------------------------------------------------------
#
# Step 1 — Start the Julia REPL:
#     Ctrl+Shift+P  (to bring up Command Palette)>  Julia: Start REPL
#
# Step 2 — Install the package and load it:
#     julia> using Pkg
#     julia> Pkg.add("CustomGaussQuadrature")
#     julia> Pkg.add("Plots")
#     julia> Pkg.add("GaussQuadrature")
#     julia> Pkg.add("SpecialFunctions")
#     julia> using CustomGaussQuadrature
#
#   If I want to, I can check the package version using:
#   julia> pkgversion(CustomGaussQuadrature)
#
# Step 3 — Run the script:
#     julia> pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
#     julia> include(joinpath(pkg_dir, "test", "test_CustomGaussQuadrature_scr.jl"))
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

# These internal functions are not part of the package's public API. The
# qualified import below makes them explicitly available in this script.
using CustomGaussQuadrature: μ_offsetvec_fn, Δ_fn, Δ′_fn, Δ_offsetvec_fn,
    Δ′_offsetvec_fn, α_offsetvec_fn, β_vec_fn, step1_fn

# Both the inbuilt package Dates and the package Printf have a format command
using Dates
date_time_now = now();
println("\n", Dates.Date(date_time_now), 
"   ", Dates.format(date_time_now, "HH:MM"))    

# pkg_dir is the package root directory, regardless of whether the
# package was loaded from a local path or from the registry.
pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
println("\n", "The package root is:")
println_wrap(pkg_dir)
print("\n")

println("This file served as a development script for")
println("CustomGaussQuadrature and was later expanded.")
println("It is a practical record of how the package was")
println("checked, with current output compared against")
println("previously saved output files.", "\n")

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
using Plots
using GaussQuadrature
using SpecialFunctions
using QuadGK
using DoubleFloats

println("------------------------------------------------------")
println("Test the moment determinants method for Step 1")
println("------------------------------------------------------", "\n")

println("------------------------------------------------------")
println("Consider which_f = ['scaled chi pdf', [0,Inf], 160]")
println("Julia's BigFloat number type has default 256 bits")
println("and corresponds to nbits <- 256 in the code for my")
println("my R package custom.gauss.quad", "\n")

println("Compare the results of low-level computations with")
println("the results of the same low-level computations in R")
println("using my R package custom.gauss.quad", "\n")

println("This comparison uses a stored weight function, so")
println("moment_fn = moment_stored_fn;")
moment_fn = moment_stored_fn;

T = BigFloat;
m = 160;
which_f = ["scaled chi pdf", [0,Inf], m]::Vector{Any};
println("T = ", T, ",    which_f = ", which_f, "\n")

n = 3;
println("n = ", n, "\n")
μ_offsetvec = μ_offsetvec_fn(BigFloat, moment_fn, which_f, n);
println("μ_offsetvec = ")
for element in μ_offsetvec
    println(element)
end
print("\n")

n = 10;
println("n = ", n, "\n")
μ_offsetvec = μ_offsetvec_fn(BigFloat, moment_fn, which_f, n);
println("Δ_fn(μ_offsetvec, n, n) = ", "\n", Δ_fn(μ_offsetvec, n, n), "\n")

println("Δ′_fn(μ_offsetvec, n, n) = ", "\n", Δ′_fn(μ_offsetvec, n, n), "\n")

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
println("Δ_offsetvec = ")
for element in Δ_offsetvec
    println(element)
end
print("\n")

Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);
println("Δ′_offsetvec = ")
for element in Δ′_offsetvec
    println(element)
end
print("\n")


α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
println("α_offsetvec = ")
for element in α_offsetvec
    println(element)
end
println("α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n) matches the results of my R code", "\n") 

β_vec = β_vec_fn(Δ_offsetvec, n)
println("β_vec = ")
for element in β_vec
    println(element)
end
println("β_vec_fn(Δ_offsetvec, n) matches the results of my R code")

println("\n","------------------------------------------------------")
println("Consider which_f = ['generalized laguerre', [0, Inf], α_gl], where α_gl=1")
println("This comparison uses a stored weight function.", "\n")


T = BigFloat;
which_f = ["generalized laguerre", [0, Inf], 1]::Vector{Any}; 
println("T = ", T, ",   which_f = ", which_f, "\n")
println("So the call laguerre_coefs(n, convert(T, 1)) uses")
println("the value α = 1 converted to type ", T, ".", "\n")

n = 10;
println("n = ", n, "\n")

μ_offsetvec = μ_offsetvec_fn(T, moment_fn, which_f, n);

Δ_offsetvec = Δ_offsetvec_fn(μ_offsetvec, n);
Δ′_offsetvec = Δ′_offsetvec_fn(μ_offsetvec, n);

α_offsetvec= α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n);
β_vec =  β_vec_fn(Δ_offsetvec, n);

println("Use the command laguerre_coefs from the GaussQuadrature package:")
println("a, b = laguerre_coefs(n, convert(T, 1))", "\n")
a, b = laguerre_coefs(n, convert(T, 1));

println("The parent function converts an Offset Array into an ordinary Array")
println("parent(α_offsetvec) - a = ")
# The parent function converts an Offset Array into an ordinary Array
for element in (parent(α_offsetvec) - a)
    println(element)
end
print("\n")

println("sqrt.(β_vec)- b[2:n] = ")
for element in (sqrt.(β_vec)- b[2:n])
    println(element)
end

println("\n","------------------------------------------------------")
println("Test of the function custom_gauss_quad_fn and print the nodes")
println("and weights in the same format as that used for the output")
println("from my R package custom.gauss.quad.")
println("This section tests Step 2 only.")
println("It checks whether the Double64 arithmetic used in Step 2")
println("works as expected by comparing the resulting nodes and")
println("weights with those obtained using BigFloat arithmetic")
println("with a 132-bit significand.", "\n")

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 33;
println("n = ", n, "\n")


println("@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

setprecision(BigFloat, 132, base=2);
println("BigFloat precision set to 132 bit significand")
println("@time nodes_BigFloat132, weights_BigFloat132 = custom_gauss_quad_fn(BigFloat, moment_fn, which_f, n);")
@time nodes_BigFloat132, weights_BigFloat132 = custom_gauss_quad_fn(BigFloat, moment_fn, which_f, n);
print("\n")

println("maximum(abs.(nodes_Double64 - nodes_BigFloat132)) = ",
    convert(Float64, maximum(abs.(nodes_Double64 - nodes_BigFloat132))))

println("maximum(abs.((nodes_Double64 - nodes_BigFloat132) ./ nodes_BigFloat132)) = ",
    convert(Float64, maximum(abs.((nodes_Double64 - nodes_BigFloat132) ./ nodes_BigFloat132))))

println("maximum(abs.((weights_Double64 - weights_BigFloat132) ./ weights_BigFloat132)) = ",
    convert(Float64, maximum(abs.((weights_Double64 - weights_BigFloat132) ./ weights_BigFloat132))), "\n")

println("The nodes and weights are printed three per row to match the")
println("presentation style used in the R package output. This makes")
println("visual comparison easier when checking saved output files.")
row_number_last = ceil(Int64, length(nodes_Double64)/3);
println("Double64 nodes:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 1]
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 2]
    @printf "%.17f \n" nodes_Double64[3*(row_number-1) + 3]
end
print("\n")
println("Double64 weights:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights_Double64[3*(row_number-1) + 1]
    @printf "%.16e " weights_Double64[3*(row_number-1) + 2]
    @printf "%.16e \n" weights_Double64[3*(row_number-1) + 3]
end

println("\n", "------------------------------------------------------")
println("Test the function custom_gauss_quad_all_fn")
println("for the scaled chi pdf.")
println("This section computes the rule using the built-in weight")
println("function interface, i.e. custom_gauss_quad_all_fn with")
println("moment_stored_fn, and then converts the nodes and weights")
println("to Float64 for printing")
println("in the same format as that used for the output from my")
println("R package custom.gauss.quad.", "\n")

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n = 33;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);
nodes = convert(Vector{Float64}, nodes_Double64);
weights = convert(Vector{Float64}, weights_Double64);
println("The nodes and weights are printed in the same format as that")
println("used for the output from my R package custom.gauss.quad.")
print("\n")
row_number_last = ceil(Int64, length(nodes)/3);

println("Double64 nodes after conversion to Float64:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes[3*(row_number-1) + 1]
    if 3*(row_number-1) + 2 ≤ n
        @printf "%.17f " nodes[3*(row_number-1) + 2]
    end
    if 3*(row_number-1) + 3 ≤ n
        @printf "%.17f \n" nodes[3*(row_number-1) + 3]
    end
end
if 3*(row_number_last-1) + 3 !== n
    println("")
end
print("\n")

println("Double64 weights after conversion to Float64:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights[3*(row_number-1) + 1]
    if 3*(row_number-1) + 2 ≤ n
        @printf "%.16e " weights[3*(row_number-1) + 2]
    end
    if 3*(row_number-1) + 3 ≤ n
        @printf "%.16e \n" weights[3*(row_number-1) + 3]
    end
end
if 3*(row_number_last-1) + 3 !== n
    println("")
end

println("\n", "------------------------------------------------------")
println("Test the function custom_gauss_quad_all_fn")
println("for the Hermite weight function.")
println("This section uses the built-in weight function interface,")
println("i.e. custom_gauss_quad_all_fn with moment_stored_fn, and")
println("compares the nodes and weights returned by this call with")
println("those returned by the")
println("hermite command from the GaussQuadrature package.")
println("BigFloat arithmetic with a 256-bit significand is used")
println("for the GaussQuadrature reference computation.")
which_f = ["hermite", [-Inf, Inf]];
println("which_f = ", which_f)

println("precision(BigFloat) = ",precision(BigFloat))
println("setprecision(BigFloat, 256, base=2) ", "\n")
setprecision(BigFloat, 256, base=2);

n= 5;
println("n = ", n, "\n")

println("Use the command hermite from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);")
@time nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);
print("\n")

println("The GaussQuadrature and CustomGaussQuadrature nodes and")
println("weights are printed below in aligned columns for direct")
println("comparison.")

println("        nodes_BigFloat (GaussQuadrature)               weights_BigFloat (GaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

println("         nodes (CustomGaussQuadrature)                weights (CustomGaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end
print("\n")

println("maximum(abs.(nodes - nodes_BigFloat)) = ",
    convert(Float64, maximum(abs.(nodes - nodes_BigFloat))))

println("maximum(abs.((weights - weights_BigFloat) ./ weights_BigFloat)) = ",
    convert(Float64, maximum(abs.((weights - weights_BigFloat) ./ weights_BigFloat))))

    
println("\n", "------------------------------------------------------")
println("Test the function custom_gauss_quad_all_fn")
println("for the chemistry example.")
println("This section uses the stored weight function interface,")
println("i.e. custom_gauss_quad_all_fn with moment_stored_fn.")
println("The resulting nodes and weights are then printed in the")
println("same format as that used in Table 2.2 of Gautschi (1983).", "\n")

which_f = ["chemistry example", [0, Inf]]::Vector{Any};
println("which_f = ", which_f, "\n")

n = 15;
println("n = ", n, "\n")
println("@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);
print("\n")

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

println("The nodes and weights are printed in the same format as that")
println("used in Table 2.2 of Gautschi (1983):")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " nodes[i]
    @printf "%.15e  \n" weights[i]
end

println("\n", "------------------------------------------------------")
println("Test the function custom_gauss_quad_all_fn")
println("for the Legendre weight function.")
println("This section uses the stored weight function interface,")
println("i.e. custom_gauss_quad_all_fn with moment_stored_fn, and")
println("compares the nodes and weights returned by this call with")
println("those returned by the")
println("legendre command from the GaussQuadrature package.")
println("BigFloat arithmetic is used for the GaussQuadrature")
println("reference computation.")
which_f = ["legendre", [-1, 1]];
println("which_f = ", which_f, "\n")

n = 5;
println("n = ", n, "\n")

println("Use the command legendre from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);")
@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);
println("precision(BigFloat) = ", precision(BigFloat))
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
print("\n")

println("The GaussQuadrature and CustomGaussQuadrature nodes and")
println("weights are printed below in aligned columns for direct")
println("comparison.")

println("        nodes_BigFloat (GaussQuadrature)               weights_BigFloat (GaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

println("         nodes (CustomGaussQuadrature)                weights (CustomGaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end
print("\n")

println("maximum(abs.(nodes - nodes_BigFloat)) = ",
    convert(Float64, maximum(abs.(nodes - nodes_BigFloat))))

println("maximum(abs.((weights - weights_BigFloat) ./ weights_BigFloat)) = ",
    convert(Float64, maximum(abs.((weights - weights_BigFloat) ./ weights_BigFloat))))

println("\n", "------------------------------------------------------")
println("Test the function custom_gauss_quad_all_fn")
println("for the generalized Laguerre weight function.")
println("This section uses the stored weight function interface,")
println("i.e. custom_gauss_quad_all_fn with moment_stored_fn, and")
println("compares the nodes and weights returned by this call with")
println("those returned by the")
println("laguerre command from the GaussQuadrature package.")
α_spec = "3.1";
which_f = ["generalized laguerre", [0, Inf], α_spec]::Vector{Any};
println("which_f = ", which_f, "\n")

α = parse(BigFloat, α_spec);
println("α = parse(BigFloat, α_spec), where α_spec = ", repr(α_spec), "\n")

n = 5;
println("n = ", n, "\n")

println("Use the command laguerre from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);")
@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
print("\n")

println("The GaussQuadrature and CustomGaussQuadrature nodes and")
println("weights are printed below in aligned columns for direct")
println("comparison.")

println("        nodes_BigFloat (GaussQuadrature)               weights_BigFloat (GaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

println("         nodes (CustomGaussQuadrature)                weights (CustomGaussQuadrature)")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes[i]
    @printf "%.33e  \n" weights[i]
end
print("\n")

println("maximum(abs.(nodes - nodes_BigFloat)) = ",
    convert(Float64, maximum(abs.(nodes - nodes_BigFloat))))

println("maximum(abs.((weights - weights_BigFloat) ./ weights_BigFloat)) = ",
    convert(Float64, maximum(abs.((weights - weights_BigFloat) ./ weights_BigFloat))))

println("\n", "------------------------------------------------------")
println("Test the function custom_gauss_quad_all_fn")
println("with upto_n = true for the scaled chi pdf.")
println("This section uses the stored weight function interface,")
println("i.e. custom_gauss_quad_all_fn with moment_stored_fn.")
println("The output consists of the Gauss quadrature rules for")
println("q = 1 up to q = n, printed one rule at a time.", "\n")

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n = 5;
println("n = ", n, "\n")

upto_n = true;
println("upto_n = ", upto_n, "\n")

println("@time nodes_upto_n, weights_upto_n = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n, upto_n);")
@time nodes_upto_n, weights_upto_n = custom_gauss_quad_all_fn(moment_stored_fn, which_f, n, upto_n);
print("\n")

println("For q = 1 the output consists of a scalar node and scalar weight.")
println("For q ≥ 2 the output consists of vectors of nodes and weights.")
print("\n")

@printf "%2d     " 1
@printf "node:   %.17f     " convert(BigFloat, nodes_upto_n[1])
@printf "weight: %.17f  \n" convert(BigFloat, weights_upto_n[1])
@printf "\n"
for q in 2:n
    @printf "%2d     \n" q
    @printf "nodes:   "
    for j in 1:q
        @printf "%.17f     " convert(BigFloat, nodes_upto_n[q][j])
    end
    @printf "\n"
    @printf "weights: "
    for j in 1:q
        @printf "%.17f    " convert(BigFloat, weights_upto_n[q][j])
    end
    @printf "\n \n"
end


println("\n", "------------------------------------------------------")
println("Try out using a new weight function, namely,")
println("the Weibull pdf with scale parameter 1")
println("and shape parameter k > 0.")
println("This section defines the corresponding moment function,")
println("computes the Gauss rule using custom_gauss_quad_all_fn,")
println("and then performs a basic validity check by comparing")
println("the empirical cdf implied by the rule with the exact")
println("Weibull cdf in a plot.", "\n")

function moment_weibull_pdf_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
@assert which_f[1] == "weibull pdf"
T_k = materialize_scalar_spec_fn(T, which_f[3])
@assert T_k > convert(T, 0)
@assert s ≥ 0
if s == 0
    return(convert(T, 1))
end
T_1 = convert(T, 1)
T_s = convert(T, s)
gamma(T_1 + (T_s/T_k))
end;

k_spec = "2.0";
k = parse(Float64, k_spec);
which_f = ["weibull pdf", [0, Inf], k_spec];
println("which_f = ", which_f, ",  k = ", k_spec, "\n")
println("k = parse(Float64, k_spec), where k_spec = ", repr(k_spec))
println("So the custom Weibull moment function uses the shape")
println("parameter value k = ", k, ".", "\n")

n = 63;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64 = ")
println("custom_gauss_quad_all_fn(moment_weibull_pdf_fn, which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_weibull_pdf_fn, which_f, n);
print("\n")
row_number_last = ceil(Int64, length(nodes_Double64)/3);
println("Double64 nodes:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 1]
    @printf "%.17f " nodes_Double64[3*(row_number-1) + 2]
    @printf "%.17f \n" nodes_Double64[3*(row_number-1) + 3]
end
print("\n")
println("Double64 weights:")
for row_number in 1:row_number_last
    @printf "[%i] " 3*(row_number-1) + 1
    @printf "%.16e " weights_Double64[3*(row_number-1) + 1]
    @printf "%.16e " weights_Double64[3*(row_number-1) + 2]
    @printf "%.16e \n" weights_Double64[3*(row_number-1) + 3]
end

# Basic validity check:
# Plot the empirical cdf corresponding to the Gauss quadrature
# rule and the exact cdf of a Weibull distribution with scale
# parameter 1 and shape parameter k.
#
# Agreement between these two curves is only a basic visual check;
# it is not a proof of accuracy of the Gauss rule.

# To include a plot pane within the VS Code editor window,
# open the command palette with Ctrl + Shift + P
# and enter the command
# Julia: Enable Plot Pane

x_vec = nodes_Double64;
prob_vec = weights_Double64;
x_lo = 0.0;
x_hi = 2.5;
p = plot_cdf_discrete_rv_fn(x_vec, prob_vec, x_lo, x_hi)
x_grid = range(x_lo, x_hi, length=200);
y_grid = weibull_cdf_fn.(x_grid, k);
plot!(p, x_grid, y_grid,
title="Weibull pdf weight fn with shape parameter k=$k & n=$n",
titlefont=font(10))
display(p)


println("\n", "------------------------------------------------------")
println("Test the Stieltjes procedure method")
println("This section compares the Stieltjes procedure with the")
println("moment-determinants method for the stored scaled chi pdf.")
println("It first compares the recurrence coefficients a_vec and")
println("b_vec, and then compares the resulting nodes and weights.", "\n")

println("Final values of a_vec and b_vec", "\n")

m = 160;
which_f = ["scaled chi pdf", [0,Inf], m];
n = 5;
println("which_f = ", which_f)
println("n = ", n, "\n")

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

println("r = ", r)

@time "step1_fn" a_vec, b_vec, μ₀, nbits = 
step1_fn(moment_stored_fn, which_f, n);

max_abs_diff_a_vec = maximum(abs.(stieltjes_a_vec - a_vec));
println("maximum(abs.(stieltjes_a_vec - a_vec)) = ", convert(Float64, max_abs_diff_a_vec))

max_abs_rel_diff_b_vec = maximum(abs.((stieltjes_b_vec - b_vec) ./ b_vec));
println("maximum(abs.((stieltjes_b_vec - b_vec) ./ b_vec)) = ", convert(Float64, max_abs_rel_diff_b_vec))

println("\n") 
println("nodes and weights", "\n")
println("The Stieltjes and moment-determinants Gauss rules are")
println("compared below using maximum absolute error for the nodes")
println("and maximum absolute relative error for the weights.", "\n")

@time "custom_gauss_quad_all_fn" nodes, weights = 
custom_gauss_quad_all_fn(moment_stored_fn, which_f, n);
print("\n")

println("maximum(abs.(stieltjes_nodes - nodes)) = ", 
convert(Float64, maximum(abs.(stieltjes_nodes - nodes))))
println("maximum(abs.((stieltjes_weights - weights) ./ weights)) = ", 
convert(Float64, maximum(abs.((stieltjes_weights - weights) ./ weights))))
print("\n")

stieltjes_nodes = convert(Vector{Float64}, stieltjes_nodes);
stieltjes_weights = convert(Vector{Float64}, stieltjes_weights);

println("           stieltjes_nodes             stieltjes_weights")
for i in 1:lastindex(stieltjes_nodes)
    @printf "%2d     " i
    @printf "%.16e     " stieltjes_nodes[i]
    @printf "%.16e  \n" stieltjes_weights[i]
end


println("\n", "------------------------------------------------------")
println("Test stieltjes_custom_gauss_quad_all_fn")
println("with upto_n = true.")
println("This section compares the Stieltjes and moment-determinants")
println("Gauss rules for q = 1 up to q = n.", "\n")

upto_n = true;
println("moment_fn = moment_stored_fn;")
moment_fn = moment_stored_fn;
println("upto_n = ", upto_n)
println("which_f = ", which_f)
println("n = ", n, "\n")

T = BigFloat;
println("T = ", T, "\n")
a = convert(T, 0);
b = Inf;

@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes_upto_n, stieltjes_weights_upto_n = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b, upto_n);  
print("\n")
@time "custom_gauss_quad_all_fn" nodes_upto_n, weights_upto_n = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, upto_n);
print("\n")

max_abs_diff_nodes = zeros(n);
for q in 1:n
    max_abs_diff_nodes[q] = maximum(abs.(stieltjes_nodes_upto_n[q] - nodes_upto_n[q]));
end
println("maximum(abs.(stieltjes_nodes_upto_n[q] - nodes_upto_n[q])) over q = 1:n = ",
    maximum(max_abs_diff_nodes))

max_abs_reldiff_weights = zeros(n);
for q in 1:n
    max_abs_reldiff_weights[q] = 
    maximum(abs.((stieltjes_weights_upto_n[q] - weights_upto_n[q]) ./ weights_upto_n[q]));
end
println("maximum(abs.((stieltjes_weights_upto_n[q] - weights_upto_n[q]) ./ weights_upto_n[q])) over q = 1:n = ",
    maximum(max_abs_reldiff_weights), "\n")

println("Convert the nodes for q = 1 up to q = n into vectors of Float64")
println("so that the node differences can be displayed clearly.", "\n")

nodes_upto_n_Float64 = Vector{Vector{Float64}}(undef, n);
nodes_upto_n_Float64[1] = [convert(Float64, nodes_upto_n[1])];
for q in 2:n
    nodes_upto_n_Float64[q] = convert(Vector{Float64}, nodes_upto_n[q]);
end

stieltjes_nodes_upto_n_Float64 = Vector{Vector{Float64}}(undef, n);
stieltjes_nodes_upto_n_Float64[1] = [convert(Float64, stieltjes_nodes_upto_n[1])];
for q in 2:n
    stieltjes_nodes_upto_n_Float64[q] = convert(Vector{Float64}, stieltjes_nodes_upto_n[q]);
end

println("stieltjes_nodes_upto_n_Float64 .- nodes_upto_n_Float64")
show(stdout, "text/plain", stieltjes_nodes_upto_n_Float64 .- nodes_upto_n_Float64)
println("\n")

println("Convert the weights for q = 1 up to q = n into vectors of Float64")
println("so that the weight differences can be displayed clearly.", "\n")

weights_upto_n_Float64 = Vector{Vector{Float64}}(undef, n);
weights_upto_n_Float64[1] = [convert(Float64, weights_upto_n[1])];
for q in 2:n
    weights_upto_n_Float64[q] = convert(Vector{Float64}, weights_upto_n[q]);
end

stieltjes_weights_upto_n_Float64 = Vector{Vector{Float64}}(undef, n);
stieltjes_weights_upto_n_Float64[1] = [convert(Float64, stieltjes_weights_upto_n[1])];
for q in 2:n
    stieltjes_weights_upto_n_Float64[q] = convert(Vector{Float64}, stieltjes_weights_upto_n[q]);
end

println("stieltjes_weights_upto_n_Float64 .- weights_upto_n_Float64")
show(stdout, "text/plain", stieltjes_weights_upto_n_Float64 .- weights_upto_n_Float64)
println("\n")


println("\n", "------------------------------------------------------")
println("Test the chemistry example weight function.")
println("This is the weight function considered by Gautschi (1983).")
println("This section compares the Stieltjes and moment-determinants")
println("Gauss rules and then prints the Stieltjes rule in the same")
println("format as Table 2.2 of Gautschi (1983).", "\n")

moment_fn = moment_stored_fn;
println("moment_fn = moment_stored_fn;")

n = 15;
which_f = ["chemistry example", [0, Inf]];
println("which_f = ", which_f)
println("number of Gauss quadrature nodes n = ", n, "\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

println("Set up the typed log-weight function and the support")
println("endpoints needed by the Stieltjes procedure.", "\n")

T = BigFloat;
println("T = ", T, "\n")
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_chemistry_fn(T, x);
μ₀ = moment_fn(T, which_f, 0);

@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b);
print("\n")

println("The Stieltjes and moment-determinants Gauss rules are")
println("compared below using maximum absolute error for the nodes")
println("and maximum absolute relative error for the weights.", "\n")

println("maximum(abs.(stieltjes_nodes - nodes)) = ",
    convert(Float64, maximum(abs.(stieltjes_nodes - nodes))))

println("maximum(abs.((stieltjes_weights - weights) ./ weights)) = ",
    convert(Float64, maximum(abs.((stieltjes_weights - weights) ./ weights))))
print("\n")

stieltjes_nodes_bf = convert(Vector{BigFloat}, stieltjes_nodes);
stieltjes_weights_bf = convert(Vector{BigFloat}, stieltjes_weights);

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to
# BigFloat, before printing using @printf.

println("Print the Stieltjes nodes and weights in the same format as")
println("Table 2.2 of Gautschi (1983):", "\n")

println("       stieltjes_nodes_bf        stieltjes_weights_bf")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " stieltjes_nodes_bf[i]
    @printf "%.15e  \n" stieltjes_weights_bf[i]
end
print("\n")

println("For easier visual comparison in saved output files, the")
println("same Stieltjes nodes and weights are also printed after")
println("conversion to Float64.", "\n")

stieltjes_nodes_Float64 = convert(Vector{Float64}, stieltjes_nodes);
stieltjes_weights_Float64 = convert(Vector{Float64}, stieltjes_weights);

println("       stieltjes_nodes_Float64   stieltjes_weights_Float64")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " stieltjes_nodes_Float64[i]
    @printf "%.15e  \n" stieltjes_weights_Float64[i]
end
print("\n")



println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Hermite
println("Hermite weight function")

setprecision(BigFloat, 256, base=2);
 
n = 15;
println("number of Gauss quadrature nodes n = ", n, "\n")

which_f = ["hermite", [-Inf, Inf]];
println("Set up the typed log-weight function and the support")
println("endpoints needed by the Stieltjes procedure.", "\n")
T = BigFloat;
println("T = ", T, "\n")
a, b = which_f[2];
lnf_typed_fn = (T, which_f, x) -> lnf_hermite_fn(T, x);
μ₀ = moment_fn(T, which_f, 0);

@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
println("r = ", r, "\n")
@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, stieltjes_a_vec, stieltjes_b_vec, a, b);  


nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
println("\n", "                nodes_BigFloat                         weights_BigFloat")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

stieltjes_nodes_bf = convert(Vector{BigFloat}, stieltjes_nodes);
stieltjes_weights_bf = convert(Vector{BigFloat}, stieltjes_weights);

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to
# BigFloat, before printing using @printf.

println("             stieltjes_nodes_bf                             stieltjes_weights_bf")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " stieltjes_nodes_bf[i]
    @printf "%.33e  \n" stieltjes_weights_bf[i]
end

print("\n")
max_abs_error_nodes = maximum(abs.(stieltjes_nodes_bf - nodes_BigFloat));
println("maximum(abs.(stieltjes_nodes_bf - nodes_BigFloat)) = ", convert(Float64, max_abs_error_nodes))

max_abs_rel_error_weights = maximum(abs.((stieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat));
println("maximum(abs.((stieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat)) = ", 
convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Generalized Laguerre
println("Generalized Laguerre weight function", "\n")

setprecision(BigFloat, 256, base=2);
n= 15;
println("number of Gauss quadrature nodes n = ", n)

T = BigFloat;
println("T = ", T, "\n")
α = convert(T, 1);
which_f = ["generalized laguerre", [0, Inf], α]::Vector{Any};
println("Set up the typed log-weight function and the support")
println("endpoints needed by the Stieltjes procedure.", "\n")
a, b = which_f[2];
println("α = convert(T, 1), so the command laguerre(n, α)")
println("uses α = 1 converted to type ", T, ".", "\n")
lnf_typed_fn = (T, which_f, x) -> lnf_laguerre_fn(T, x, which_f[3]);

nodes_BigFloat, weights_BigFloat = laguerre(n, α);

μ₀ = moment_fn(T, which_f, 0);
@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_typed_fn, which_f, a, b);
println("r = ", r, "\n")

@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀,stieltjes_a_vec, stieltjes_b_vec, a, b);

print("\n")
println("                nodes_BigFloat                         weights_BigFloat")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " nodes_BigFloat[i]
    @printf "%.33e  \n" weights_BigFloat[i]
end
print("\n")

stieltjes_nodes_bf = convert(Vector{BigFloat}, stieltjes_nodes);
stieltjes_weights_bf = convert(Vector{BigFloat}, stieltjes_weights);

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to
# BigFloat, before printing using @printf.

println("             stieltjes_nodes_bf                             stieltjes_weights_bf")
for i in 1:n
    @printf "%2d   " i
    @printf "%.33e   " stieltjes_nodes_bf[i]
    @printf "%.33e  \n" stieltjes_weights_bf[i]
end
print("\n")

max_abs_error_nodes = maximum(abs.(stieltjes_nodes_bf - nodes_BigFloat));
println("maximum(abs.(stieltjes_nodes_bf - nodes_BigFloat)) = ", convert(Float64, max_abs_error_nodes))

max_abs_rel_error_weights = maximum(abs.((stieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat));
println("maximum(abs.((stieltjes_weights_bf - weights_BigFloat) ./ weights_BigFloat)) = ", 
convert(Float64, max_abs_rel_error_weights))

println("\n", "NOTE: VS Code diff intra-line highlighting may be unreliable for the last few changed lines.")

