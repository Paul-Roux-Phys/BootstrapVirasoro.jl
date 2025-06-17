#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

#====
This file presents the whole user-interface of the package.
====#

module BootstrapVirasoro

using EllipticFunctions, Printf

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
export Spectrum,
    ChannelSpectrum,
    ChannelSpectra,
    add!,
    remove!,
    fields,
    hasdiagonals

# Bootstrap equations and solver.
export BootstrapSystem,
    evaluate_blocks!,
    compute_linear_system!,
    solve!,
    write_csv

# print complex numbers with 5 digits
function format_complex(z::Complex{<:Real})
    real_str = @sprintf("%.5e", real(z))
    imag_str = @sprintf("%.5e", abs(imag(z)))
    sign = imag(z) < 0 ? "-" : "+"
    buf = real(z) > 0 ? " " : ""
    return "$buf$real_str $sign $(imag_str)im"
end

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
include("CFTData/CentralCharges.jl")
include("CFTData/ConformalDimensions.jl")
include("CFTData/Fields.jl")

# Conformal blocks
# AbstractBlocks serve as an interface to all types of blocks.
include("ConformalBlocks/AbstractBlocks.jl")

# Linear bootstrap equations
include("BootstrapEquations/Spectrum.jl")
include("BootstrapEquations/StructureConstants.jl")
include("BootstrapEquations/LinearSystem.jl")

end
