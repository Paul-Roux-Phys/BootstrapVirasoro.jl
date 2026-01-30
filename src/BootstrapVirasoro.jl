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
export LeftRight, Field, swap_lr, total_dimension, spin
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
export evaluate_correlation

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
import Printf: Format, format

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
struct Channels{T}
    dict::Dict{Symbol, T}
end
const Chans = Channels

Base.getindex(x::LeftRight, s::Symbol) = getfield(x, s)
Base.getindex(s::Channels, ch::Symbol) = s.dict[ch]
function Base.setindex!(s::Channels, value, ch::Symbol)
    s.dict[ch] = value
end
function Base.setproperty!(s::Channels, ch::Symbol, value)
    if ch in (:s, :t, :u)
        s.dict[ch] = value
    else
        Base.setfield!(s, value, ch)
    end
end
Channels(s::T, t, u) where {T} = Channels{T}(Dict(:s => s, :t => t, :u => u))
Channels{T}(s, t, u) where {T} = Channels{T}(Dict{Symbol,T}(:s => s, :t => t, :u => u))
Channels(t::NTuple{3,T}) where {T} = Channels(t[1], t[2], t[3])
Channels(s) = Channels(s, s, s)
Channels(f::Function) = Channels(f(:s), f(:t), f(:u))
Channels(v::Dict) = Channels{v.parameters[2]}(v)
Base.length(c::Channels) = length(c.dict)
function Base.getproperty(c::Channels, s::Symbol)
    s in (:s, :t, :u) && return getfield(c, :dict)[s]
    return getfield(c, s)
end
Channels(::Type{T}) where {T} = Channels{T}(Dict{Symbol, T}())
Base.keys(c::Channels) = Base.keys(c.dict)
Base.haskey(c::Channels, k::Symbol) = Base.haskey(c.dict, k)

@static if VERSION >= v"1.9"
    max_thread_id() = Threads.maxthreadid()
else
    max_thread_id() = Threads.nthreads()
end

# Include generic interface for in-place arithmetics on bigfloats
include("MutableArithmetics_ComplexBigFloat.jl")

# Implementations
include("1_CFTData.jl")
include("2_Correlation.jl")
include("3_Block.jl")
include("4_bootstrapsystem.jl")
include("LoopModels.jl")

end # module BootstrapVirasoro
