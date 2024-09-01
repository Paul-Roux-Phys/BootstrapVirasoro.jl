#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module JuliVirBootstrap

using Latexify, # print outputs in latex format
    Base.Threads

# The exported methods and types are found in the following included files.
# These files also document the exported methods
include("JuliVirBootstrap/CFTData.jl")
include("JuliVirBootstrap/ConformalBlocks.jl")

end
