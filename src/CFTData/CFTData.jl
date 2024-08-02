#===========================================================================================

CFTData.jl provides types representing
central charges and fields in 2D CFTs with Virasoro symmetry.

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module CFTData

using Latexify;
import ..JuliVirBootstrap: left, right

export CentralCharge,
       ConformalDimension,
       Field, spin, swap_lr

include("central_charges.jl")
include("conformal_dimensions.jl")
include("fields.jl")

end # end module
