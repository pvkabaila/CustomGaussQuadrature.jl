#  conversion_to_Double64_review_scr.jl
# This script is for investigation, not for automated testing.

using Pkg
Pkg.add("DoubleFloats")
Pkg.add("Printf")

using DoubleFloats
using Printf

mu0 = convert(Double64, 1)
@printf "%.32e     " mu0

mu0 = convert(Double64, 1.0)
@printf "%.128e     " mu0

Dbl_mu0 = convert(Double64, 1.1)
@printf "%.128e     " mu0
BigFloat_mu0 = convert(BigFloat, 1.1)
@printf "%.16e     " Dbl_mu0 - BigFloat_mu0

mu0 = 1.1
parse(BigFloat, string(1.1))

# Exact integer inputs such as mu0 = 1 are safe; this conversion stays exact.
mu0 = convert(Double64, 1)
@printf "%.128e     " mu0

k = 1.1
Dbl_k = parse(Double64, string(1.1))
BigFloat_k = parse(BigFloat, string(k))
Dbl_k - BigFloat_k