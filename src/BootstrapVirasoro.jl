#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

#====
This file presents the whole user-interface of the package.
====#

module BootstrapVirasoro

using EllipticFunctions, BarnesDoubleGamma, DataFrames, CSV

using GenericLinearAlgebra: qr

# CFT Data: Central charges, conformal dimensions, fields.
export LeftRight,
    CentralCharge,
    CC,
    ConformalDimension,
    CD,
    P_rs,
    δrs,
    Field,
    spin,
    isdiagonal,
    isdegenerate,
    swap_lr,
    shift,
    reflect,
    total_dimension

# Conformal blocks and correlations.
export Correlation,
    Block,
    qfromx,
    xfromq,
    qfromτ,
    evaluate_series,
    evaluate_series_der,
    evaluate,
    evaluate_der

# Spectra
export Spectrum, ChannelSpectrum, ChannelSpectra, add!, remove!, fields, hasdiagonals

# Bootstrap equations and solver.
export BootstrapSystem, evaluate_blocks!, compute_linear_system!, solve!, write_csv

"""
        LeftRight{T}
Left and right pairs of objects. Can be accessed with
`obj[:left]` and `obj[:right]`.
"""
const LeftRight{T} = Tuple{T,T} # left and right pairs of things

function Base.getindex(x::LeftRight, s::Symbol)
    s === :left && return x[1]
    s === :right && return x[2]
    error("tried to access pair element other than 1, 2, :left or :right")
end

# The exported methods and types are found in the following included files.
# The files also document the exported methods
# CFT Data: central charges, conformal dimensions, fields
include("CFTData/central_charges.jl")
include("CFTData/conformal_dimensions.jl")
include("CFTData/fields.jl")

# Conformal blocks
# AbstractBlocks serve as an interface to all types of blocks.
include("ConformalBlocks/abstract_blocks.jl")

# Linear bootstrap equations
include("BootstrapEquations/Spectrum.jl")
include("BootstrapEquations/structure_constants.jl")
include("BootstrapEquations/linear_system.jl")

end
