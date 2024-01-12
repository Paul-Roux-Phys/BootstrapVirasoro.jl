#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram 
Nivesvivat

===========================================================================================#

module JuliVirBootstrap

#===========================================================================================
Central charges and fields
===========================================================================================#
include("CFTData.jl")
using .CFTData
export CentralCharge, Field


#===========================================================================================
Conformal blocks
===========================================================================================#
include("ConformalBlocks.jl")
using .FourPointBlocksSphere, .OnePointBlocksTorus
#export 


#===========================================================================================
Special functions
===========================================================================================#
include("SpecialFunctions.jl")
export log_double_gamma, double_gamma


end