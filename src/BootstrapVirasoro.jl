#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

#====
This file presents the API (user-interface) of the package.
====#

module BootstrapVirasoro

using EllipticFunctions, BarnesDoubleGamma, DataFrames, CSV, PrettyTables

using Printf: @sprintf
using GenericLinearAlgebra: qr

#======================================================================

CFT Data: Central charges, conformal dimensions, fields.

=======================================================================#
export CentralCharge, CC
export ConformalDimension, CD
export LeftRight, Field, swap_lr, shift, reflect, total_dimension

struct CentralCharge{T}
    β::T
    B::T
    b::T
    c::T
    n::T
end
# convenience alias
const CC = CentralCharge

struct ConformalDimension{T}
    c::CentralCharge{T}
    P::T
    p::T
    δ::T
    Δ::T
    r::Union{Rational,Int}   # r must be rational
    s::Union{T,Rational,Int} # s can be an arbitrary number to represent a
    isKac::Bool              # diagonal field. Still needs to be exact in
    degenerate::Bool         # case of a non-diagonal field
end
# convenience alias
const CD = ConformalDimension

struct LeftRight{T}
    left::T
    right::T
end
const LR = LeftRight
Base.getindex(x::LeftRight, s::Symbol) = getfield(x, s)

struct Field{T}
    c::CC{T}
    dims::LeftRight{ConformalDimension{T}}
    r::Union{Rational,Int}
    s::Union{T,Rational,Int}
    diagonal::Bool
    degenerate::Bool
    isKac::Bool
end

#=======================================================
Conformal blocks and correlations.
========================================================#
export Correlation, Corr, Co
export Block, IBlock

abstract type Correlation{T} end
const Corr = Correlation # type alias
const Co = Correlation # type alias

abstract type Block{T} end # general conformal block. Can be interchiral, non-chiral or chiral

# # Spectra
# export Spectrum, ChannelSpectrum, ChannelSpectra, add!, remove!, fields, hasdiagonals

# # Bootstrap equations and solver.
# export BootstrapSystem, evaluate_blocks!, compute_linear_system!, solve!

# The implementations are found in the following included files.
# CFT Data: central charges, conformal dimensions, fields
include("CFTData.jl")
include("ConformalBlocks/correlations.jl")


# Conformal blocks
# AbstractBlocks serve as an interface to all types of blocks.
include("ConformalBlocks/abstract_blocks.jl")

# Linear bootstrap equations
include("BootstrapEquations/Spectrum.jl")
# include("BootstrapEquations/structure_constants.jl")
# include("BootstrapEquations/linear_system.jl")

end

"""
    CentralCharge(param = value)

Type representing a central charge, with precision given by the type of `value`.
The supported `param` are `c`, `β`, `b`, `B`.

# Examples

```jldoctest
julia> c = CentralCharge(β=big"0.1"+big"0.2"*im); c.β ≈ big"0.1" + big"0.2"*im
true

julia> c = CentralCharge(c = 0.7);

julia> c.b ≈ -0.0 + 0.894427190999916im
true

julia> c.n ≈ 1.6180339887498953 + 0.0im
true
```
"""

"""
    ConformalDimension(c; param=value, r=0, s=0)
    
Type representing a conformal dimension, with precision given by the central charge.
The supported parameters are `Δ`, `δ`, `P`, `p`, or the Kac indices `r, s`.

# Examples

```jldoctest
julia> c = CentralCharge(β = 0.3im);

julia> d1 = ConformalDimension(c, δ=0.5); d1.Δ ≈ 3.800277777777777 + 0.0im
true

julia> d1.P ≈ 0.7071067811865476
true

julia> d2 = ConformalDimension(c, r=3//2, s=2//3)
ConformalDimension{ComplexF64} with Kac indices r = 3//2, s = 2//3
```
"""

"""
        LeftRight{T}
Left and right pairs of objects. Can be accessed with
`obj[:left]` and `obj[:right]`.
"""

"""
    Field(c, parameter=value; r, s, diagonal=false)

Type for representing a non-chiral field.

# keyword arguments:

- `r::Rational`,`s::Rational`. By convention ``V_{(r,s)}`` has
left and right momenta ``(P_{(r,s)}, P_{(r,-s)})``.
- `diagonal::Bool`: set to `true` to get a diagonal field;

# Examples

```jldoctest
julia> c = CentralCharge(β = big"0.5");

julia> V = Field(c, r=0, s=1)
Diagonal Field{Complex{BigFloat}}, dim = Δ_{0, 1}

julia> (V.dims[:left].Δ, V.dims[:right].Δ)
(0.4375 + 0.0im, 0.4375 + 0.0im)

julia> V.dims[:left].P ≈ 1
true

julia> V.dims[:right].p ≈ -im
true

julia> V2 = Field(c, :P, 0.42, diagonal=true); V2.diagonal
true
```
"""

"""
    Correlation(args...; Δmax=10.)

Abstract type for holding data relevant to the computation of a correlation function:
* external fields
* Coefficients ``R_{m,n}``, possibly left and right and for different channels
* Coefficients ``R^{\\text{reg}}_{m, n}`` when ``R_{m, n}`` vanishes
* Coefficients ``C^N_{m, n}``, possibly left and right and for different channels
up to `N=Δmax`.

It is also possible to access the central charge, left or right parts of the correlation if it is non chiral, etc. See examples.

# Examples

```jldoctest
julia> c = CentralCharge(β = sqrt(2));

julia> V1 = Field(c, r=2, s=3//2);

julia> co = Correlation(V1, V1, V1, V1, 10)
NonChiralCorrelation{ComplexF64,NTuple{4, Field{ComplexF64}}} with external fields
< V_{(2, 3//2)} V_{(2, 3//2)} V_{(2, 3//2)} V_{(2, 3//2)} >

julia> co._Rmn[:left][:s][2, 4] ≈ 8.37053571428573e-7 - 0.0im
true

julia> co.c.c ≈ -2
true

julia> co.fields
(V_{(2, 3//2)}, V_{(2, 3//2)}, V_{(2, 3//2)}, V_{(2, 3//2)})

julia> co[:left]
ChiralCorrelation{ComplexF64} with external dimensions
(Δ_{2, 3//2}, Δ_{2, 3//2}, Δ_{2, 3//2}, Δ_{2, 3//2})
```
"""

"""
        Block(co, chan, d or V; Δmax=10., interchiral=false, der=false)

Object representing a block in a channel. Can be a chiral, derivative, non-chiral, or interchiral block, depending on the types of the input correlation and channel dimension `d` or channel field `V`. Can be evaluated at a position with [`evaluate`](@ref).

# Example

```jldoctest
julia> c = CentralCharge(β = sqrt(2));

julia> V1 = Field(c, r=2, s=3//2);

julia> co = Correlation(V1, V1, V1, V1, 10);

julia> V = Field(c, r=2, s=1//2); b = Block(co, :s, V)
Non chiral factorised block for the correlation
< V_{(2, 3//2)} V_{(2, 3//2)} V_{(2, 3//2)} V_{(2, 3//2)} >
channel: s, V_{(2, 1//2)}

julia> b = Block(co[:left], :s, V.dims[:left])
Chiral block for the correlation
< Δ_{2, 3//2} Δ_{2, 3//2} Δ_{2, 3//2} Δ_{2, 3//2} >
Channel: s, Δ_{2, 1//2}
```
"""
