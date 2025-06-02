#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

#====
This file presents the whole user-interface of the package.
====#

module BootstrapVirasoro

using EllipticFunctions,
    Printf

using GenericLinearAlgebra: qr


# CFT Data: Central charges, conformal dimensions, fields.
export LeftRight,
    CentralCharge, CC,
    ConformalDimension, CD, P_rs,
    Field, spin, isdiagonal, isdegenerate,
    swap_lr, shift,
    total_dimension

# Conformal blocks and correlations.
export Correlation,
    Block,
    qfromx, xfromq,
    qfromτ,
    evaluate_series, evaluate_series_der,
    evaluate, evaluate_der

# Spectra
export Spectrum,
    ChannelSpectrum,
    ChannelSpectra,
    add!, remove!

# Bootstrap equations and solver.
export BootstrapSystem,
    compute_linear_system!,
    solve!

# The exported methods and types are found in the following included files.
# These files also document the exported methods
include("CFTData.jl")
include("ConformalBlocks.jl")
include("spectra.jl")
include("bootstrap_equations.jl")
include("linear_solver.jl")

end
