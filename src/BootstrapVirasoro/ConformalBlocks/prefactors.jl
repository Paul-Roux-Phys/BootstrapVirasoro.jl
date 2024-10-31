include("prefactors/four_point.jl")
include("prefactors/one_point.jl")

blockprefactor_chiral(b::BlockChiral, x) = block_prefactor_chiral(b.corr.dims, b, x)