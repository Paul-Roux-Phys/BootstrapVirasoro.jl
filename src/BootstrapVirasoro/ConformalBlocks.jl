using Memoization,
    EllipticFunctions
import BarnesDoubleGamma: digamma_reg

export Correlation,
    ExtFields,
    FourFields,
    OneField,
    Block,
    BlockChiral,
    BlockNonChiral,
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

abstract type Block{T} end

Block(c, chan, V::Field, Nmax::Int) = BlockNonChiral(c, chan, V, Nmax)
Block(c, chan, V::Field, lr::Symbol, Nmax; der=false) = 
    BlockChiral(c, chan, V.dims[lr], lr, Nmax, der=der)
Block(c, chan, d::ConformalDimension, lr, Nmax; der=false) =
    BlockChiral(c, chan, d, lr, Nmax, der=der)
Block(c, chan, V::Field, lr::Symbol) = Block(c, chan, V, lr, c.Nmax)
Block(c, chan, V::Field) = Block(c, chan, V, c.Nmax)

include("ConformalBlocks/residues.jl")
include("ConformalBlocks/Correlations.jl")
include("ConformalBlocks/BlocksSeries.jl")
include("ConformalBlocks/prefactors.jl")
include("ConformalBlocks/logarithmic.jl")
include("ConformalBlocks/evaluate.jl")
