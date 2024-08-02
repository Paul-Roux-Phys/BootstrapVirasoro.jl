#===========================================================================================

ConformalBlocks.jl contains modules that compute Virasoro four-point conformal blocks on the
sphere and Virasoro one-point conformal blocks on the torus.

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#


"""
Computation of conformal blocks via Zamolodchikov recursion.
"""
module ConformalBlocks

using Memoization,
    EllipticFunctions
using ..CFTData
import ..JuliVirBootstrap: left, right

export FourPointCorrelation, 
    OnePointCorrelation,
    FourPointBlock,
    OnePointBlock,
    evaluate_series,
    evaluate_chiral

abstract type Correlation{T} end
abstract type Block{T} end

include("Correlations/Correlations.jl")
include("BlocksSeries/BlocksSeries.jl")
include("prefactors/prefactors.jl")
include("evaluate/evaluate.jl")

end # end module
