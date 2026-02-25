
using CustomGaussQuadrature
using Printf
using Plots
using GaussQuadrature
using SpecialFunctions



m = 160;
which_f = ["scaled chi pdf", [0,Inf], m];
n = 33;


#  @__DIR__ gives the folder of the file currently being executed
#  .. = parent folder (package root),
#  then into src,
#  then to stjieltjes_lnf_stored_scr.jl
include(joinpath(@__DIR__, "..", "src", "stjieltjes_lnf_stored_scr.jl"));

#-----------------------------------------------------------------
# Compare nodes and weight obtained using the Stjieltjes procedure
# with those obtained using moment determinants

nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

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



which_f = ["Hermite", [-Inf,Inf]];
n = 10;
#  @__DIR__ gives the folder of the file currently being executed
#  .. = parent folder (package root),
#  then into src,
#  then to stjieltjes_lnf_stored_scr.jl
include(joinpath(@__DIR__, "..", "src", "stjieltjes_lnf_stored_scr.jl"));
#-----------------------------------------------------------------
# Compare nodes and weight obtained using the Stjieltjes procedure
# with those obtained using moment determinants

nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n);

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