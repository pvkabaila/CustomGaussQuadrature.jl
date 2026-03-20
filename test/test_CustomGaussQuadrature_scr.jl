# This script tests the functions for computing the Gauss
# rule using the (a) moment determinant method and (b) 
# Stieltjes procedure.
# The scripts under test reside in
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
#     julia> include("test/test_CustomGaussQuadrature_scr.jl")
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
# Step 3 — Run the script:
#     julia> pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
#     julia> include(joinpath(pkg_dir, "test", "test_CustomGaussQuadrature_scr.jl"))
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
#     julia> Pkg.add("Plots")
#     julia> Pkg.add("GaussQuadrature")
#     julia> Pkg.add("SpecialFunctions")
#     julia> using CustomGaussQuadrature
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

# pkg_dir is the package root directory, regardless of whether the
# package was loaded from a local path or from the registry.
pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
println("\n", "The package root is:")
println_wrap(pkg_dir)
print("\n")

println("The simplest way to compare the output from this test,")
println("stored as a plain text file, with a previous output")
println("from this test is as follows. Open both plain text files")
println("in VS Code. Then on the OPEN EDITORS menu on the upper left")
println("right click on the previous output > Select for Compare.")
println("Then right click on the text file for the latest output >")
println("Compare with Selected. Both plain text files are compared")
println("in the same window, with corresponding line numbers from")
println("both plain text files.")
println("\n", "\n")

println("2026 3 17 I used Claude Opus 4.6 to check the safety of installing")
println("the EXTENSION Persistent Highlighter and it reported that no") 
println("dangerous, suspicious or questionable code was found.") 
println("To always highlight with the color Light Yellow use") 
println("Ctrl + Shift + P to bring up the Command Palette and") 
println("then use Persistent Highlighter: Add Custom Color Highlight")
println("and always pick Light Yellow.", "\n")

using Printf
using Plots
using GaussQuadrature
using SpecialFunctions
using QuadGK
using DoubleFloats

#****************************************************
# Test of gauss_quad_moment_dets_scr.jl
#****************************************************
# Try which_f = ["scaled chi pdf", [0,Inf], m]
# Julia's BigFloat number type corresponds to 
# 256 bits i.e.
# nbits <- 256
# in my R code.

println("------------------------------------------------------")
println("Test the moment determinants method")
println("------------------------------------------------------", "\n")

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
print("\n")
# α_offsetvec_fn(Δ_offsetvec, Δ′_offsetvec, n) matches the results of my R code  ✔

β_vec = β_vec_fn(Δ_offsetvec, n)
println("β_vec = ")
for element in β_vec
    println(element)
end
# β_vec_fn(Δ_offsetvec, n) matches the results of my R code  ✔

println("\n", "------------------------------------------------------")
#----------------------------------------------------
# Try which_f = ["Generalized Laguerre", [0, Inf], α_GGL], where α_GGL=1

T = BigFloat;
which_f = ["Generalized Laguerre", [0, Inf], 1]::Vector{Any}; 
println("T = ", T, ",   which_f = ", which_f, "\n")

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

println("\n", "------------------------------------------------------")
#--------------------------------------------------------------
# Test of the function custom_gauss_quad_fn and print the nodes
# and weights in the same format as that used for the output
# from my R package custom.gauss.quad

println("Test of the function custom_gauss_quad_fn and print the nodes")
println("and weights in the same format as that used for the output")
println("from my R package custom.gauss.quad.", "\n")

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

abs_error_nodes = convert(Vector{Float64}, abs.(nodes_Double64 - nodes_BigFloat132));
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", max_abs_error_nodes)

abs_rel_error_nodes = convert(Vector{Float64}, (abs.((nodes_Double64 - nodes_BigFloat132) ./ nodes_BigFloat132)));
max_abs_rel_error_nodes = maximum(abs_rel_error_nodes);
println("max_abs_rel_error_nodes = ", max_abs_rel_error_nodes)

abs_rel_error_weights = convert(Vector{Float64}, (abs.((weights_Double64 - weights_BigFloat132) ./ weights_BigFloat132)));
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights, "\n")

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
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# scaled chi pdf
which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n = 33;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{Float64}, nodes_Double64);
weights = convert(Vector{Float64}, weights_Double64);
# Print the nodes and weights in the same format as that used for the output
# from my R package custom.gauss.quad
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
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Hermite
which_f = ["Hermite", [-Inf, Inf]];
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

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

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

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# chemistry example
which_f = ["chemistry example", [0, Inf]]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 15;
println("n = ", n, "\n")
println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
print("\n")

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);

#---------------------------------------------------------------------------
# Print the nodes and weights in the same format as that used in Table 2.2
# of Gautschi (1983)
println("Print the nodes and weights in the same format as that used in")
println("Table 2.2 of Gautschi (1983):")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " nodes[i]
    @printf "%.15e  \n" weights[i]
end

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Legendre
which_f = ["Legendre", [-1, 1]];
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("Use the command legendre from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);")
@time nodes_BigFloat, weights_BigFloat = legendre(BigFloat, n);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
println("precision(BigFloat) = ", precision(BigFloat))
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
print("\n")

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

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function custom_gauss_quad_all_fn
# Generalized Laguerre
α = 3.1;
which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
println("which_f = ", which_f)

println("α = ", α, "\n")

n= 5;
println("n = ", n, "\n")

println("Use the command laguerre from the GaussQuadrature package:")
println("@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);")
@time nodes_BigFloat, weights_BigFloat = laguerre(n, α);
print("\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);
nodes = convert(Vector{BigFloat}, nodes);
weights = convert(Vector{BigFloat}, weights);
print("\n")

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

error_nodes = nodes - nodes_BigFloat;
abs_error_nodes = abs.(error_nodes);
max_abs_error_nodes = maximum(abs_error_nodes);
println("max_abs_error_nodes = ", convert(Float64, max_abs_error_nodes))

rel_error_weights = (weights - weights_BigFloat) ./ weights_BigFloat;
abs_rel_error_weights = abs.(rel_error_weights);
max_abs_rel_error_weights = maximum(abs_rel_error_weights);
println("max_abs_rel_error_weights = ", convert(Float64, max_abs_rel_error_weights))

println("\n", "------------------------------------------------------")
#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=false and extra_check=true

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights =")
println("custom_gauss_quad_all_fn(which_f, n, false, true);")
@time nodes_Double64, weights_Double64, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, false, true);
print("\n")

println("nodes_Double64 = ")
for element in nodes_Double64
    println(element)
end
print("\n")

println("weights_Double64 = ")
for element in weights_Double64
    println(element)
end
print("\n")

println("max_abs_error_nodes = ", max_abs_error_nodes)
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights)

println("\n", "------------------------------------------------------")
#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=false

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n, true);")
@time nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n, true);
print("\n")


@printf "%2d     " 1
    @printf "node:   %.17f     " convert(BigFloat, nodes[1])
    @printf "weight: %.17f  \n" convert(BigFloat, weights[1])
    @printf "\n"
for i in 2:n
    @printf "%2d     \n" i
    @printf "nodes:   "
    for j in 1:i
        @printf "%.17f     " convert(BigFloat, nodes[i][j])
    end
    @printf "\n"
    @printf "weights: "
    for j in 1:i
        @printf "%.17f    " convert(BigFloat, weights[i][j])
    end
    @printf "\n \n"
end

println("\n", "------------------------------------------------------")
#------------------------------------------------------------
# Test the function custom_gauss_quad_all_fn with 
# upto_n=true and extra_check=true

println("Test the function custom_gauss_quad_all_fn with") 
println("upto_n=true and extra_check=true", "\n")

which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
println("which_f = ", which_f, "\n")

n= 5;
println("n = ", n, "\n")

println("@time nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights =")
println("custom_gauss_quad_all_fn(which_f, n, true, true);")
@time nodes, weights, max_abs_error_nodes, max_abs_rel_error_weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, true, true);
print("\n")


@printf "%2d     " 1
    @printf "node:   %.17f     " convert(BigFloat, nodes[1])
    @printf "weight: %.17f  \n" convert(BigFloat, weights[1])
    @printf "\n"
for i in 2:n
    @printf "%2d     \n" i
    @printf "nodes:   "
    for j in 1:i
        @printf "%.17f     " convert(BigFloat, nodes[i][j])
    end
    @printf "\n"
    @printf "weights: "
    for j in 1:i
        @printf "%.17f    " convert(BigFloat, weights[i][j])
    end
    @printf "\n \n"
end

println("max_abs_error_nodes = ",  max_abs_error_nodes)
println("max_abs_rel_error_weights = ", max_abs_rel_error_weights)


###################################################
# Use, for example,
# moment_fn = moment_weibull_pdf_fn
###################################################
println("\n", "------------------------------------------------------")
println("Try out using a new weight function, namely,") 
println("the Weibull pdf, with scale parameter 1")
println("and shape parameter γ > 0, weight function")

function moment_weibull_pdf_fn(::Type{T}, which_f, r::Integer) where {T<:AbstractFloat}
@assert which_f[1] == "weibull pdf"
γ = which_f[3]
@assert γ > 0
T_γ = parse(T, string(γ))
@assert r ≥ 0
if r == 0
    return(convert(T, 1))
end
T_1 = convert(T, 1)
T_r = convert(T, r)
gamma(T_1 + (T_r/T_γ))
end;

println("moment_fn = moment_weibull_pdf_fn;")
moment_fn = moment_weibull_pdf_fn;

γ = 2.0;
which_f = ["weibull pdf", [0, Inf], γ];
println("which_f = ", which_f, ",  γ=", γ, "\n")

n = 63;
println("n = ", n, "\n")

println("@time nodes_Double64, weights_Double64 = ")
println("custom_gauss_quad_all_fn(moment_fn, which_f, n);")
@time nodes_Double64, weights_Double64 = custom_gauss_quad_all_fn(moment_fn, which_f, n);
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

#-------------------------------------------------
# Plot the cdf corresponding to the Gauss quadrature 
# rule and the actual cdf of a Weibull distribution 
# with scale parameter λ = 1 and shape parameter γ 
# (k > 0)

# To include a plot pane within the VS code editor window,
# open the command palette with Ctrl + Shift + P
# and enter the command
# Julia: Enable Plot Pane

x_vec = nodes_Double64;
prob_vec = weights_Double64;
x_lo = 0.0;
x_hi = 2.5;
p = plot_cdf_discrete_rv_fn(x_vec, prob_vec, x_lo, x_hi) 
x_grid = range(x_lo, x_hi, length=200);
y_grid = weibull_cdf_fn.(x_grid, γ);
plot!(x_grid, y_grid, 
title="Weibull pdf weight function with scale parameter λ=1.0 and γ=$γ",
 titlefont=font(10))
display(p)

println("\n", "------------------------------------------------------")
println("Test the Stieltjes procedure method")
println("------------------------------------------------------", "\n")
#****************************************************
#  Test of gauss_quad_stieltjes_scr.jl
#***************************************************
println("Test gauss_quad_stieltjes_scr.jl", "\n") 
# Final values of a_vec and b_vec
println("Final values of a_vec and b_vec", "\n") 

println("moment_fn = moment_stored_fn;")
moment_fn = moment_stored_fn;

m = 160;
lnf_fn = x -> ln_scaled_chi_pdf_fn(T, x, m);
println("scaled chi pdf weight function, with m = ", m)
T = BigFloat;
a = convert(T,0);
b = Inf;
which_f = ["scaled chi pdf", [0,Inf], m];
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);
n = 5;
println("number of Gauss quadrature nodes n = ", n)

@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r)

@time "step1_fn" a_vec, b_vec, μ₀, nbits = 
step1_fn(moment_fn, which_f, n);

max_abs_diff_a_vec = maximum(abs.(stieltjes_a_vec - a_vec));
println("maximum(abs.(stieltjes_a_vec - a_vec)) = ", convert(Float64, max_abs_diff_a_vec))

max_abs_rel_diff_b_vec = maximum(abs.((stieltjes_b_vec - b_vec) ./ b_vec));
println("maximum(abs.((stieltjes_b_vec - b_vec) ./ b_vec)) = ", convert(Float64, max_abs_rel_diff_b_vec))

println("\n") 
# nodes and weights
println("nodes and weights", "\n") 

@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);

@time "custom_gauss_quad_all_fn" nodes, weights = 
custom_gauss_quad_all_fn(moment_fn, which_f, n);
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
#-----------------------------------------
# Test stieltjes_custom_gauss_quad_all_fn
# with upto_n = true
#-----------------------------------------
println("Test stieltjes_custom_gauss_quad_all_fn")
println("with upto_n = true", "\n")

upto_n = true;
println("upto_n = ", upto_n, "\n")

T = BigFloat;
a = convert(T,0);
b = Inf;

@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes_upto_n, stieltjes_weights_upto_n = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b, upto_n);
print("\n")
@time "custom_gauss_quad_all_fn" nodes_upto_n, weights_upto_n = 
custom_gauss_quad_all_fn(moment_fn, which_f, n, upto_n);
print("\n")

max_abs_diff_nodes = zeros(n);
for q in 1:n
    max_abs_diff_nodes[q] = maximum(abs.(stieltjes_nodes_upto_n[q] - nodes_upto_n[q]));
end
println("maximum(max_abs_diff_nodes) = ", maximum(max_abs_diff_nodes))

max_abs_reldiff_weights = zeros(n);
for q in 1:n
    max_abs_reldiff_weights[q] = 
    maximum(abs.((stieltjes_nodes_upto_n[q] - nodes_upto_n[q]) ./ nodes_upto_n[q]));
end
println("maximum(max_abs_reldiff_weights) = ", maximum(max_abs_reldiff_weights), "\n")

#---------------------------------------------------------
# Conversion of nodes_upto_n into a vector of Float64 vectors

nodes_upto_n_Float64 = Vector{Vector}(undef, n);
nodes_upto_n_Float64[1] = [convert(Float64, nodes_upto_n[1])];
for q in 2:n
    nodes_upto_n_Float64[q] = convert(Vector{Float64}, nodes_upto_n[q]);
end
# println("nodes_upto_n_Float64")
# nodes_upto_n_Float64

#---------------------------------------------------------
# Conversion of stieltjes_nodes_upto_n into a vector of Float64 vectors

stieltjes_nodes_upto_n_Float64 = Vector{Vector}(undef, n);
stieltjes_nodes_upto_n_Float64[1] = [convert(Float64, stieltjes_nodes_upto_n[1])];
for q in 2:n
    stieltjes_nodes_upto_n_Float64[q] = convert(Vector{Float64}, stieltjes_nodes_upto_n[q]);
end
# println("stieltjes_nodes_upto_n_Float64")
# stieltjes_nodes_upto_n_Float64

println("stieltjes_nodes_upto_n_Float64 .- nodes_upto_n_Float64")
show(stdout, "text/plain", stieltjes_nodes_upto_n_Float64 .- nodes_upto_n_Float64)
println("\n")

#---------------------------------------------------------
# Conversion of weights_upto_n into a vector of Float64 vectors

weights_upto_n_Float64 = Vector{Vector}(undef, n);
weights_upto_n_Float64[1] = [convert(Float64, weights_upto_n[1])];
for q in 2:n
    weights_upto_n_Float64[q] = convert(Vector{Float64}, weights_upto_n[q]);
end
# println("weights_upto_n_Float64")
# weights_upto_n_Float64


#---------------------------------------------------------
# Conversion of stieltjes_weights_upto_n into a vector of Float64 vectors

stieltjes_weights_upto_n_Float64 = Vector{Vector}(undef, n);
stieltjes_weights_upto_n_Float64[1] = [convert(Float64, stieltjes_weights_upto_n[1])];
for q in 2:n
    stieltjes_weights_upto_n_Float64[q] = convert(Vector{Float64}, stieltjes_weights_upto_n[q]);
end
# println("stieltjes_weights_upto_n_Float64")
# stieltjes_weights_upto_n_Float64

println("stieltjes_weights_upto_n_Float64 .- weights_upto_n_Float64")
show(stdout, "text/plain", stieltjes_weights_upto_n_Float64 .- weights_upto_n_Float64)


println("\n", "------------------------------------------------------")
#--------------------------------------------------
# The weight function considered by Gautschi (1983)
println("The weight function considered by Gautschi (1983).")
println("In other words, the chemistry example weight function.", "\n")

n = 15;
println("number of Gauss quadrature nodes n = ", n)

which_f = ["chemistry example", [0, Inf]];
nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

lnf_fn = x -> lnf_chemistry_fn(T, x);
T = BigFloat;
a = convert(T,0);
b = Inf;

which_f = ["chemistry example", [0, Inf]];
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1);

print("\n")
@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);
print("\n")

# Beware!
# If @printf is applied to a Double64 number then it
# first converts it to a Float64 number before
# printing, irrespective of the number of digits specified.
# A way to fix this is to first convert the Double64 to 
# BigFloat, before printing using @printf.
stieltjes_nodes_bf = convert(Vector{BigFloat}, stieltjes_nodes);
stieltjes_weights_bf = convert(Vector{BigFloat}, stieltjes_weights);

#---------------------------------------------------------------------------
# Print the nodes and weights in the same format as that used in Table 2.2
# of Gautschi (1983)
println("Print the nodes and weights in the same format as that used in")
println("Table 2.2 of Gautschi (1983):", "\n")

println("       stieltjes_nodes_bf        stieltjes_weights_bf")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " stieltjes_nodes_bf[i]
    @printf "%.15e  \n" stieltjes_weights_bf[i]
end
print("\n")

#---------
stieltjes_nodes_Float64 = convert(Vector{Float64}, stieltjes_nodes);
stieltjes_weights_Float64 = convert(Vector{Float64}, stieltjes_weights);
#---------------------------------------------------------------------------
# Print the nodes and weights in the same format as that used in Table 2.2
# of Gautschi (1983)
println("       stieltjes_nodes_Float64   stieltjes_weights_Float64")
for i in 1:n
    @printf "%2d     " i
    @printf "%.15e     " stieltjes_nodes_Float64[i]
    @printf "%.15e  \n" stieltjes_weights_Float64[i]
end

#-----------------------------------------------



println("\n", "------------------------------------------------------")
#------------------------------------------------
# Test the function stieltjes_custom_gauss_quad_all_fn
# Hermite
println("Hermite weight function")

setprecision(BigFloat, 256, base=2);
 
n = 15;
println("number of Gauss quadrature nodes n = ", n, "\n")

which_f = ["Hermite", [-Inf, Inf]];
lnf_fn = lnf_hermite_fn;
a = -Inf;
b = Inf;
T = BigFloat;
μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1); 

@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r, "\n")
@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);


nodes_BigFloat, weights_BigFloat = hermite(BigFloat, n);
println("\n", "                nodes_BigFloat                         weights_BigFloat")
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
stieltjes_nodes_bf = convert(Vector{BigFloat}, stieltjes_nodes);
stieltjes_weights_bf = convert(Vector{BigFloat}, stieltjes_weights);
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
α = convert(T, 1);
which_f = ["Generalized Laguerre", [0, Inf], α]::Vector{Any};
a = convert(T,0);
b = Inf;
println("Generalized Laguerre, with α = ", α, "\n")
lnf_fn = x -> lnf_laguerre_fn(x, α);

nodes_BigFloat, weights_BigFloat = laguerre(n, α);

μ₀, μ₁ = μ_offsetvec_fn(T, moment_fn, which_f, 1); 
@time "stieltjes_a_vec_b_vec_final_fn" stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
println("r = ", r, "\n")

@time "stieltjes_custom_gauss_quad_all_fn" stieltjes_nodes, stieltjes_weights = 
stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);

print("\n")
println("                nodes_BigFloat                         weights_BigFloat")
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
stieltjes_nodes_bf = convert(Vector{BigFloat}, stieltjes_nodes);
stieltjes_weights_bf = convert(Vector{BigFloat}, stieltjes_weights);
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
