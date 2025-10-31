#===========================================================================================

Written by Paul Roux, adapting and extending a Python code written by Sylvain Ribault &
Rongvoram Nivesvivat.

===========================================================================================#

module BootstrapVirasoro

#=========================================================================
API
=========================================================================#
export CentralCharge, CC
export ConformalDimension, CD
export LeftRight, Field, swap_lr, reflect, total_dimension, spin
export Channels, @channels, Chans
export Correlation, Corr, Co, Correlation4, Correlation1
export ChiralCorrelation, ChiralCorrelation4, ChiralCorrelation1
export NonChiralCorrelation, NonChiralCorrelation4, NonChiralCorrelation1
export ChiralBlock, CBlock, NonChiralBlock, NCBlock
export LinearCombinationBlock, LCBlock, Block
export ChannelSpectrum, ChanSpec, add!, remove!
export StructureConstants, StrCst, SC, find_normalised
export BootstrapSystem, evaluate_blocks!, compute_linear_system!, solve!
export solve_bootstrap, solve_bootstrap_manyP

#=========================================================================
Third-party functions
=========================================================================#
import EllipticFunctions: jtheta2, jtheta3, ellipticK, etaDedekind
import BarnesDoubleGamma: DoubleGamma, gamma, digamma_reg
import DataFrames: DataFrame, groupby, select, Not
import CSV
using LaTeXStrings
import PrettyTables: pretty_table, @crayon_str, Highlighter, Tables.columntable
import Random
import Printf: @sprintf

#========================================================================
Package-wide definitions
========================================================================#
struct LeftRight{T}
    left::T
    right::T
end
const LR = LeftRight
const NB_CHANNELS = 3
const CHANNELS = (:s, :t, :u)
mutable struct Channels{T}
    s::T
    t::T
    u::T
end
const Chans = Channels

Base.getindex(x::LeftRight, s::Symbol) = getfield(x, s)

Base.getindex(s::Channels, ch::Symbol) = getfield(s, ch)
Base.setindex!(s::Channels, value, ch::Symbol) = setfield!(s, ch, value)
Channels(s::T, t, u) where {T} = Channels{T}(s, t, u)
Channels(t::NTuple{3,T}) where {T} = Channels{T}(t[1], t[2], t[3])
Channels(s) = Channels(s, s, s)
Channels(f::Function) = Channels(f(:s), f(:t), f(:u))
Base.length(c::Channels) = 3

# Implementations
include("1_CFTData.jl")
include("2_Correlation.jl")
include("3_Block.jl")
include("4_bootstrapsystem.jl")
include("LoopModels.jl")

end # module BootstrapVirasoro
