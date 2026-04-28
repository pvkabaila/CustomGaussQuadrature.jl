using CustomGaussQuadrature
using Test

@testset "Driver-owned T in Stieltjes user path" begin

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))

@eval Main begin
    which_f = ["Hermite", [-Inf, Inf]]
    n = 3
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))

stored_nodes = copy(Main.stieltjes_nodes)
stored_weights = copy(Main.stieltjes_weights)
stored_a_vec = copy(Main.stieltjes_a_vec)
stored_b_vec = copy(Main.stieltjes_b_vec)
stored_nbits = Main.stieltjes_nbits
stored_r = Main.r

@eval Main begin
    function lnf_user_hermite_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
        @assert which_f[1] == "Hermite"
        T_x = convert(T, x)
        -(T_x * T_x)
    end

    which_f = ["Hermite", [-Inf, Inf]]
    n = 3
    lnf_typed_fn = lnf_user_hermite_fn
    mu0 = sqrt(big(pi))
end
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

@test all(isapprox.(Main.stieltjes_nodes, stored_nodes; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(Main.stieltjes_weights, stored_weights; rtol=1.0e-24, atol=1.0e-28))
@test all(isapprox.(Main.stieltjes_a_vec, stored_a_vec; rtol=big"1e-24", atol=big"1e-70"))
@test all(isapprox.(Main.stieltjes_b_vec, stored_b_vec; rtol=big"1e-24", atol=big"1e-70"))
@test Main.stieltjes_nbits == stored_nbits
@test Main.stieltjes_r == stored_r

end