"""
Functionality relevant for bootstraping loop models.
"""
module LoopModels

export InterchiralBlock, IBlock, shift

using ..BootstrapVirasoro
import ..BootstrapVirasoro: getfields, total_dimension, reflect
import ..BootstrapVirasoro: reflect
import SpecialFunctions: gamma
import BarnesDoubleGamma: DoubleGamma

struct InterchiralBlock{T} <: LinearCombinationBlock{T}
    chan_field::Field{T}
    blocks::Vector{Block{T}}
    coeffs::Vector{T}
    chan::Symbol
end
const IBlock = InterchiralBlock

# not the same as V.r, V.s because it returns (0, -2βP) also for degenerate fields
# the latter returns integer (r, s) for the degenerate fields.
function indices(V)
    if V.degenerate
        0, -2*V.c.β*V[:left].P
    else
        V.r, V.s
    end
end


""" s → s+shift"""
shift(d::CD, shift) = CD(d.c, r = d.r, s = d.s + shift)

"""
        shift(V, i)

Shift the field:
- s -> s+i if !V.diagonal
"""
function shift(V::Field, i)
    if V.diagonal
        Field(shift(V.dims[:left], i))
    else
        Field(V.c, r=V.r, s=V.s + i)
    end
end


""""
        shift_C(V1, V2, V3)

Ratio ``C_{123}/C_{123^++}``.
See arXiv:2411.17262 (4.13a).
"""
function shift_C123(V1, V2, V3)
    β = V1.c.β
    r1, s1 = indices(V1)
    r2, s2 = indices(V2)
    r3, s3 = indices(shift(V3, 1))

    res = (-1)^(Int(2r2 + max(2r1, 2r2, 2r3, r1 + r2 + r3)))
    res *= (β + 0im)^(-4 / β^2 * s3)
    res *= prod(
        gamma(
            1 // 2 +
            1 // 2 * abs(pm1 * r1 + pm2 * r2 - r3) +
            1 / (2 * β^2) * (pm1 * s1 + pm2 * s2 - s3),
        ) / gamma(
            1 // 2 +
            1 // 2 * abs(pm1 * r1 + pm2 * r2 + r3) +
            1 / (2 * β^2) * (pm1 * s1 + pm2 * s2 + s3),
        ) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

"""
        shift_C122(V1, V2)

Ratio ``C_{122}/C_{12^++2^++}``.
"""
function shift_C122(V1, V2)
    β = V1.c.β
    r1, s1 = indices(V1)
    r2, s2 = indices(shift(V2, 1))

    res = (-1)^(Int(2r2))
    res *= β^(-8 / β^2 * s2)
    res *= prod(
        gamma(
            1 // 2 + 1 // 2 * abs(r1 + 2pm2 * r2) -
            1 / (β^2 * 2) * (pm2 * s1 + 2s2 + pm1),
        ) / gamma(
            1 // 2 +
            1 // 2 * abs(r1 + 2pm2 * r2) +
            1 / (β^2 * 2) * (pm2 * s1 + 2s2 + pm1),
        ) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

""""
        shift_B(V)

Ratio ``B_{V}/B_{V^++}``.
See arXiv:2411.17262 (4.13b).
"""
function shift_B(V)
    β = V.c.β
    r, s = indices(shift(V, 1))
    res = (-1)^(2r)
    res *= β^(-8inv(β)^2 * s)
    res *= prod(
        gamma(a + r - s / β^2) / gamma(a + r + s / β^2) for
        a in (0, 1, inv(β^2), 1 - inv(β^2))
    )
end

"""ratio ``D_{V}/D_{V^++}`` for four-point structure constants"""
function shift_D(Vs::NTuple{4,Field}, V)
    shift_C123(Vs[1], Vs[2], V) * shift_C123(Vs[3], Vs[4], V) / shift_B(V)
end

""" ratio ``D_{V}/D_{V^++}`` for torus one-point structure constants"""
function shift_D(Vs::NTuple{1,Field}, V)
    shift_C122(Vs[1], V) / shift_B(V)
end

function Base.length(b::IBlock)
    length(b.blocks)
end

function my_mod2(s)
    s = s % 2
    return s > 1 ? s - 2 : s
end

function _InterchiralBlock(co::Correlation{T}, chan, V, Δmax, _shift) where {T}
    Δmax === missing && (Δmax = co.Δmax)
    @assert real(V.c.β^2) > 0 "interchiral blocks are implemented for Re(β^2) > 0"
    blocks = []
    shifts = []
    extfields = getfields(co, chan)
    V0 = V
    if !V.diagonal
        V0 = Field(V.c, r = V.r, s = my_mod2(V.s))
    end
    Vshift = V0
    D = 1
    while real(total_dimension(Vshift)) <= real(Δmax) && spin(Vshift) < real(Δmax)
        push!(blocks, NCBlock(co, chan, Vshift, Δmax))
        push!(shifts, D)
        D /= shift_D(extfields, Vshift)
        Vshift = shift(Vshift, _shift)
    end
    if !(V.r isa Int && V.s isa Int)
        Vshift = shift(V0, -_shift)
        D = shift_D(extfields, Vshift)
        while real(total_dimension(Vshift)) <= real(Δmax) && spin(Vshift) < real(Δmax)
            push!(blocks, NCBlock(co, chan, Vshift, Δmax))
            push!(shifts, D)
            Vshift = shift(Vshift, -_shift)
            D *= shift_D(extfields, Vshift)
        end
    end

    InterchiralBlock{T}(V, blocks, shifts, chan)
end

function InterchiralBlock(
    co::Co,
    chan,
    V::Field,
    Δmax = missing,
    shift = 2;
    parity = 0,
)
    if parity == 0 || V.diagonal || V.s == 0 || V.s == 1
        _InterchiralBlock(co, chan, V, Δmax, shift)
    else
        _InterchiralBlock(co, chan, V, Δmax, shift) +
        parity * _InterchiralBlock(co, chan, reflect(V), Δmax, shift)
    end
end

function InterchiralBlock(
    co::Correlation1,
    V::Field,
    Δmax = missing,
    shift = 2;
    parity = 0,
)
    InterchiralBlock(co, :τ, V, Δmax, shift, parity = parity)
end

reflect(b::IBlock) =
    IBlock(b.blocks[1].corr, b.chan, reflect(b.chan_field), b.blocks[1].Δmax)

function Base.show(io::IO, ::MIME"text/plain", b::IBlock)
    svals = [B.chan_field.s for B in b.blocks]
    V = b.chan_field
    ns = real.(round.(-V.s .+ svals)) ./ 2
    shift_range = Int(minimum(ns)):1:Int(maximum(ns))
    if isempty(b.blocks)
        print(io, "Empty interchiral block")
    elseif !V.diagonal || V.degenerate
        println(io, "Interchiral block with channel $(b.chan) and fields")
        print(
            io,
            "{ V_{$(b.chan_field.r), $(b.chan_field.s) + 2s}, s ∈ $(shift_range) }",
        )
    else
        println(io, "Interchiral block with channel $(b.chan) and fields")
        print(io, "V_{P = $(V[:left].P) + n/β}, n ∈ $(shift_range)}")
    end
end

function Base.show(io::IO, b::IBlock)
    svals = [B.chan_field.s for B in b.blocks]
    V = b.blocks[1].chan_field
    ns = real(round.(-V.s .+ svals)) ./ 2
    shift_range = Int(minimum(ns)):1:Int(maximum(ns))
    if isempty(b.blocks)
        print(io, "Empty interchiral block")
    elseif !b.blocks[1].chan_field.diagonal || b.blocks[1].chan_field.degenerate
        print(io, "G^($(b.chan))({ ")
        print(
            io,
            "V_{$(b.chan_field.r), $(b.chan_field.s) + 2s}, s ∈ $(shift_range) })",
        )
    else
        print(io, "G^($(b.chan))({ ")
        print(io, "V_{P = $(V.dims[:left].P) + n/β}, n ∈ $(shift_range)})")
    end
end

function Cref(V₁, V₂, V₃, DG)
    β = V₁.c.β
    r₁, s₁ = indices(V₁)
    r₂, s₂ = indices(V₂)
    r₃, s₃ = indices(V₃)

    return prod(
        1 / DG(
            (β + 1 / β) / 2 +
            β / 2 * abs(pm₁ * r₁ + pm₂ * r₂ + pm₃ * r₃) +
            1 / 2 / β * (pm₁ * s₁ + pm₂ * s₂ + pm₃ * s₃),
        ) for pm₁ in (-1, 1) for pm₂ in (-1, 1) for pm₃ in (-1, 1)
    )
end

Cref(V₁, V₂, V₃) = Cref(V₁, V₂, V₃, DoubleGamma(V₁.c.β))

function Bref(DG, c, r, s, reg = 1 / big(10^15))
    β = c.β
    if r % 1 == 0 && s % 1 == 0
        s += reg
    end
    π = oftype(c.β, Base.π) # π in the correct precision
    return (-1)^(round(Int, r * s)) / 2 / sin(π * (r % 1 + s)) /
           sin(π * (r + s / β^2)) / prod(
        DG(β + pm1 * β * r + pm2 * s / β) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

function Bref(DG, c, P)
    β = c.β
    prod(1 / DG(β^pm1 + pm2 * 2P) for pm1 in (-1, 1) for pm2 in (-1, 1))
end

function Bref(V::Field, DG, reg = 1 / big"10"^15)
    c = V.c
    if V.diagonal
        return Bref(DG, c, V.dims[:left].P)
    else
        return Bref(DG, c, V.r, V.s, reg)
    end
end

function compute_reference(co::Correlation4, V::Field, chan, DG)
    V₁, V₂, V₃, V₄ = getfields(co, chan)
    return Cref(V₁, V₂, V, DG) * Cref(V₃, V₄, V, DG) / Bref(V, DG)
end

function compute_reference(co::Correlation1, V::Field, chan, DG)
    V₁ = getfields(co, chan)
    return Cref(V₁, V, V, DG) / Bref(V, DG)
end

compute_reference(b::Block, DG) = compute_reference(b.corr, b.chan_field, b.chan, DG)

end # end BootstrapVirasoro.LoopModels
