import BarnesDoubleGamma: digamma_reg

using EllipticFunctions
using Memoization

export Correlation,
    ExtFields,
    FourFields,
    OneField,
    Block,
    BlockChiral,
    BlockNonChiral,
    GeneralizedBlock,
    qfromx,
    xfromq,
    evaluate_series,
    evaluate

const ExtDimensions{T} = Tuple{Vararg{ConformalDimension{T}}} # tuples of dimensions
const FourDimensions{T} = NTuple{4, ConformalDimension{T}}
const OneDimension{T} = Tuple{ConformalDimension{T}}

const ExtFields{T} = Tuple{Vararg{Field{T}}} # tuples of fields
const FourFields{T} = NTuple{4, Field{T}}
const OneField{T} = Tuple{Field{T}}

"""
    Correlation{T}

Abstract type for holding data relevant to the computation of a correlation function:
* external fields
* Coefficients ``R_{m,n}``, possibly left and right and for different channels
* Coefficients ``R^{\\text{reg}}_{m, n}`` when ``R_{m, n}`` vanishes
* Coefficients ``C^N_{m, n}``, possibly left and right and for different channels

It is also possible to access the central charge, left or right parts of the correlation if
it is non chiral, etc. See examples.

# Examples

```jldoctest
julia> c = CentralCharge(β = sqrt(2))
c = -2.0000000000000027 + 0.0im, β = -1.4142135623730951 - 0.0im

julia> V1 = Field(c, r=2, s=3//2)
Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2

julia> co = Correlation(V1, V1, V1, V1, 10)
CorrelationNonChiral with external fields
Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2
Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2
Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2
Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2


julia> co._Rmn[:left][:s]
Dict{Tuple{Int64, Int64}, ComplexF64} with 17 entries:
  (1, 2)  => -0.046875-0.0im
  (2, 5)  => -0.0-0.0im
  (1, 4)  => 5.93736e12+0.0im
  (3, 2)  => -0.0-0.0im
  (4, 1)  => 0.000145397+0.0im
  (2, 1)  => 0.143555-0.0im
  (10, 1) => 9.57818e-12+0.0im
  (4, 2)  => -0.0-0.0im
  (2, 2)  => -5.93736e12-0.0im
  (1, 6)  => 1.78429e-20+0.0im
  (2, 3)  => -0.0-0.0im
  (8, 1)  => 2.23115e-9+0.0im
  (2, 4)  => 8.37054e-7-0.0im
  (5, 2)  => -0.0-0.0im
  (6, 1)  => 5.3357e-7+0.0im
  (1, 8)  => 2.84603e-23+0.0im
  (1, 10) => 7.64062e-26+0.0im

julia> co.c
c = -2.0000000000000027 + 0.0im, β = -1.4142135623730951 - 0.0im

julia> co.fields
(Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2, Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2, Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2, Non-diagonal Field{ComplexF64}
left: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
right: ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = -3//2)

julia> co[:left]
Chiral correlation function with external dimensions
ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2
ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2

```
"""
abstract type Correlation{T} end

abstract type GeneralizedBlock{T} end # generalised block, can be interchiral or conformal
abstract type Block{T} <: GeneralizedBlock{T} end # conformal block
abstract type BlockInterchiral{T} <: GeneralizedBlock{T} end # interchiral block (chiral or not)

include("ConformalBlocks/residues.jl")
include("ConformalBlocks/Correlations.jl")
include("ConformalBlocks/BlocksSeries.jl")
include("ConformalBlocks/prefactors.jl")
include("ConformalBlocks/logarithmic.jl")
include("ConformalBlocks/evaluate.jl")
include("ConformalBlocks/InterchiralBlocks.jl")

Block(c::Correlation, chan::Symbol, V::Field, Nmax::Int) =
    BlockNonChiral(c, chan, V, Nmax)
Block(c::Correlation, chan, V::Field, lr::Symbol, Nmax; der=false) = 
    BlockChiral(c, chan, V.dims[lr], lr, Nmax, der=der)
Block(c, chan, d::ConformalDimension, lr, Nmax; der=false) =
    BlockChiral(c, chan, d, lr, Nmax, der=der)
Block(c, chan, V::Field, lr::Symbol) = Block(c, chan, V, lr, c.Nmax)
Block(c, chan, d::ConformalDimension) = BlockChiral(c, chan, d, c.Nmax)
Block(c, chan, V::Field) = Block(c, chan, V, c.Nmax)

Nmax(Δmax::ConformalDimension, V::Field) =
    ceil(Int, Δmax.Δ - minimum(abs(V.dims[lr].Δ) for lr in (:left, :right)))

Block(c, chan, V::Field, Δmax::ConformalDimension) = Block(c, chan, V, Nmax(Δmax, V))
Block(c, chan, V::Field, lr::Symbol, Δmax::ConformalDimension; der=false) =
    Block(c, chan, V, Nmax(Δmax, V))
Block(c, chan, d::ConformalDimension, lr, Δmax::ConformalDimension; der=false) =
    BlockChiral(c, chan, d, lr, Nmax(Δmax, V), der=der)

function GeneralizedBlock(c, chan, V::Field, Δmax::ConformalDimension; interchiral=true)
    if interchiral
        BlockInterchiralNonChiral(c, chan, V, Δmax)
    else
        Block(c, chan, V, Δmax)
    end
end
