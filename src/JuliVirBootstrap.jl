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
export FourPointBlockSphere, F_four_point_sphere

using .OnePointBlocksTorus
export OnePointBlockTorus, F_one_point_torus


#===========================================================================================
Special functions
===========================================================================================#
include("SpecialFunctions.jl")
export log_double_gamma, double_gamma


end
