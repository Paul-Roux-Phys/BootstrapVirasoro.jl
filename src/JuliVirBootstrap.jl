#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module JuliVirBootstrap

#===========================================================================================
Special functions
===========================================================================================#
include("SpecialFunctions.jl")
using .SpecialFunctions
export Barnes_G
export log_double_Gamma, double_Gamma

#===========================================================================================
Central charges and fields
===========================================================================================#
include("CFTData.jl")
using .CFTData
export CentralCharge
export Field

#===========================================================================================
Correlation functions
===========================================================================================#
include("CorrelationFunctions.jl")
using .FourPointCorrelationFunctions
export FourPointCorrelation

using .OnePointCorrelationFunctions
export OnePointCorrelation

#===========================================================================================
Conformal blocks
===========================================================================================#
include("ConformalBlocks.jl")
using .FourPointBlocksSphere
export FourPointBlockSphere
export block_chiral, block_non_chiral

using .OnePointBlocksTorus
export OnePointBlockTorus

end
