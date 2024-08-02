#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module JuliVirBootstrap

using Latexify, # print outputs in latex format
    Base.Threads

const left  = 1
const right = 2

#===========================================================================================
Central charges and fields
===========================================================================================#
include("CFTData/CFTData.jl")
using .CFTData
export CentralCharge,
       ConformalDimension,
       Field, spin, swap_lr

#===========================================================================================
Conformal blocks
===========================================================================================#
include("ConformalBlocks/ConformalBlocks.jl")
using .ConformalBlocks
export FourPointCorrelation, 
    OnePointCorrelation,
    FourPointBlock,
    OnePointBlock,
    evaluate_series,
    evaluate_chiral
# export FourPointBlockSphere
# export block_chiral, block_non_chiral

# using .OnePointBlocksTorus
# export OnePointBlockTorus

end
