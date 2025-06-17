"""
    Spectrum(fields, Δmax; interchiral=false)

Abstract type for representing a CFT spectrum.

# examples

```jldoctest
julia> c = CentralCharge(β=1/(big"0.8"+big"0.1"*im));

julia> V1 = Field(c, r=2, s=0);

julia> co = Correlation(V1, Δmax=10.);

julia> fields = [Field(c, r=r, s=s) for r in 2:30 for s in -1+1//r:1//r:1 if r*s% 1 == 0];

julia> s = Spectrum(fields, 10.0)
Non-diagonal:
(2, -1//2), (2, 0), (2, 1//2), (2, 1)
(3, -2//3), (3, -1//3), (3, 0), (3, 1//3), (3, 2//3), (3, 1)

julia> add!(s, Field(c, r=1, s=0)); s
Non-diagonal:
(1, 0)
(2, -1//2), (2, 0), (2, 1//2), (2, 1)
(3, -2//3), (3, -1//3), (3, 0), (3, 1//3), (3, 2//3), (3, 1)
```
"""
abstract type Spectrum{T} end

include("BootstrapEquations/Spectrum.jl")

"""
TODO
"""
struct StructureConstants{T}
    constants::Channels{Dict{Field{T},T}}
    errors::Channels{Dict{Field{T},T}}
end

include("BootstrapEquations/StructureConstants.jl")

struct BootstrapMatrix{T}
    unknowns::Channels{Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
end

"""
TODO
"""
mutable struct BootstrapSystem{T,U<:ChannelSpectrum{T}}
    positions::Vector{T}
    positions_cache::Channels{Vector{LRPositionCache{T}}}        # positions at which eqs are evaluated
    spectra::Channels{U}    # channel spectra
    block_values::Channels{Dict{Field{T}, Vector{T}}} # all blocks evaluated at all positions
    matrix::BootstrapMatrix{T}  # matrix of equations
    consts::StructureConstants{T}
end

include("BootstrapEquations/LinearSystem.jl")
