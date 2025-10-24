"""
        Block(co, chan, V, Δmax=co.Δmax) # 4pt 𝕊^2 block
        Block(co, V, Δmax=co.Δmax)       # 1pt 𝕋^2 block

Given a 
"""
abstract type Block{T} end # general conformal block. Can be interchiral, non-chiral or chiral

struct ChiralBlock{T} <: Block{T}
    corr::CCo{T}
    chan_dim::CD{T}
    coeffs::Vector{T}
    coeffs_der::Vector{T}
    missing_terms::Vector{T}
    chan::Symbol
    Δmax::Int
end
const CBlock = ChiralBlock

abstract type NonChiralBlock{T} <: Block{T} end
const NCBlock = NonChiralBlock

struct FactorizedBlock{T} <: NonChiralBlock{T}
    chan_field::Field{T}
    corr::Corr
    cblocks::LR{CBlock{T}}
    Δmax::Int
    chan::Symbol
end

struct LogBlock{T} <: NonChiralBlock{T}
    cblocks::LR{CBlock{T}} # blocks for V_(r, s) x V_(r, -s)
    cblocks_op::LR{CBlock{T}} # blocks for V_(r, -s) x V_(r, s)
    cblocks_der::Union{LR{CBlock{T}},Nothing}
    corr::Corr
    chan_field::Field{T}
    R::T
    Rbar::T
    ell::T
    Δmax::Int
    nbzeros::Int
    chan::Symbol
end

abstract type LinearCombinationBlock{T} <: Block{T} end
const LCBlock = LinearCombinationBlock

struct GenericLCBlock{T} <: LCBlock{T}
    chan_field::Field{T}
    blocks::Vector{Block{T}}
    coeffs::Vector{T}
    chan::Symbol
end

#=============================================================================
Chiral Blocks
=================================++++========================================#
function series_H(d::CD{T}, Δmax, CNmn) where {T}
    P = d.P
    P2 = P^2
    isKac = d.isKac
    r, s = d.r, d.s

    coeffs = Matrix{T}(undef, Δmax + 2, Δmax + 2)
    all_mns = union([CNmn.keys[N] for N = 1:Δmax]...)
    for (m, n) in all_mns
        if isKac && m == r && n == s
            coeffs[m, n] = -inv(4CNmn.δs[r, s])
        else
            coeffs[m, n] = inv(P2 - CNmn.δs[m, n])
        end
    end

    H = zeros(T, Δmax + 1)
    H[1] = one(T)
    for N = 1:Δmax
        for (m, n) in CNmn.keys[N]
            (N > Δmax || m > Δmax || n > Δmax) && continue
            H[N+1] += CNmn[N, m, n] * coeffs[m, n]
        end
    end

    return H
end

function series_H_der(d::CD{T}, Δmax, CNmn) where {T}
    P = d.P
    P2 = P^2

    H = zeros(T, Δmax + 1)
    coeffs = Dict(
        (m, n) => -2P / (P2 - CNmn.δs[m, n])^2 for
        (m, n) in union([CNmn.keys[N] for N = 1:Δmax]...)
    )

    for N = 1:Δmax
        for (m, n) in CNmn.keys[N]
            (N > Δmax || m > Δmax || n > Δmax) && continue
            H[N+1] += CNmn[N, m, n] * coeffs[(m, n)]
        end
    end

    return H
end

function ChiralBlock(co::CCo{T}, chan, d, Δmax = missing, der = false) where {T}
    Δmax === missing && (Δmax = co.Δmax)
    CNmn = getCNmn(co, chan)
    coeffs = series_H(d, Δmax, CNmn)
    coeffs_der = Vector{T}()
    der && (coeffs_der = series_H_der(d, Δmax, CNmn))
    if !d.degenerate
        missing_terms = Vector{T}()
    else
        r, s = d.r, d.s
        missing_terms = [
            (N > 0 && (r, s) in CNmn.keys[N]) ? CNmn[N, r, s] : zero(T) for N = 0:Δmax
        ]
    end

    CBlock{T}(co, d, coeffs, coeffs_der, missing_terms, chan, Δmax)
end

function ChiralBlock(co::ChiralCorrelation1, d::CD = CD(), Δmax=missing, der=false)
    CBlock(co, :τ, d, Δmax, der)
end

getRmn(b::CBlock) = getRmn(b.corr, b.chan)
getRmnreg(b::CBlock) = getRmnreg(b.corr, b.chan)
getCNmn(b::CBlock) = getCNmn(b.corr, b.chan)
getc(b::CBlock) = b.corr.c

function Base.show(io::IO, ::MIME"text/plain", b::CBlock)
    println(io, "Chiral block for the correlation")
    println(io, b.corr)
    println(io, "Channel: $(b.chan), $(b.chan_dim)")
end

function Base.show(io::IO, b::CBlock)
    println(io, "F^($(b.chan))($(b.chan_dim))")
end

#================================================================================
Non Chiral blocks
================================================================================#
function FactorizedBlock(co::NCCo{T}, chan, V, Δmax) where {T}
    bl = CBlock(co[:left], chan, V[:left], Δmax)
    br = CBlock(co[:right], chan, V[:right], Δmax)
    FactorizedBlock{T}(V, co, LR(bl, br), Δmax, chan)
end

function islogarithmic(V::Field, _::Correlation)
    r, s = V.r, V.s
    !V.diagonal && r * s != 0 && (r % 1 == s % 1 == 0) && r > 0 && s > 0
end

islogarithmic(b::NCBlock) = islogarithmic(b.chan_field, b.corr)

function isaccidentallynonlogarithmic(co::NCCo, chan, V)
    !islogarithmic(V, co) && return false
    return Rmn_zero_order(V.r, V.s, getfields(co[:left], chan)) > 0
end

isaccidentallynonlogarithmic(b::NCBlock) =
    isaccidentallynonlogarithmic(b.corr, b.chan, b.chan_field)

function LogBlock(co::NCCo{T}, chan, V, Δmax) where {T}
    V_op = swap_lr(V)
    VV = V.s > 0 ? (V, V_op) : (V_op, V) # (V_(r, s > 0), V_(r, -s))
    r, s = VV[1].r, VV[1].s
    left, left_op, right, right_op = Tuple(
        CBlock(co[lr], chan, v.dims[lr], Δmax, false) for lr in (:left, :right)
        for v in VV
    )
    R, Rbar, l = zero(T), zero(T), zero(T)
    cl, cr = co[:left], co[:right]
    if isaccidentallynonlogarithmic(co, chan, V)
        chiral_blocks_der = nothing
        if V.s > 0
            R = getRmnreg(cl, chan)[r, s]
            Rbar = getRmnreg(cr, chan)[r, s]
        end
        nbzeros = Rmn_zero_order(r, s, left.corr.fields)
        l = zero(T)
    else
        leftder = CBlock(cl, chan, VV[2][:left], Δmax, true)
        rightder = CBlock(cr, chan, VV[1].dims.right, Δmax, true)
        chiral_blocks_der = LR(leftder, rightder)
        if V.s > 0
            R = getRmn(cl, chan)[r, s]
            Rbar = getRmn(cr, chan)[r, s]
            l = ell(getfields(co, chan), r, VV[1].s)
        end
        nbzeros = 0
    end
    LogBlock{T}(
        LR(left, right),
        LR(left_op, right_op),
        chiral_blocks_der,
        co,
        V,
        R,
        Rbar,
        l,
        Δmax,
        nbzeros,
        chan,
    )
end

reflect(b::NonChiralBlock) =
    FactorizedBlock(b.corr, b.channel, reflect(b.chan_field), b.Δmax)

function NonChiralBlock(co::NCCo, chan, V, Δmax=missing)
    Δmax === missing && (Δmax = co.Δmax)
    if islogarithmic(V, co)
        LogBlock(co, chan, V, Δmax)
    else
        FactorizedBlock(co, chan, V, Δmax)
    end
end
NonChiralBlock(co::NonChiralCorrelation1, V, Δmax=missing) = NCBlock(co, :τ, V, Δmax)

Base.getindex(b::NCBlock, s::Symbol) = b.cblocks[s]

function ell(V::NTuple{4,Field}, r, s)
    β = V[1].c.β
    βsq = β^2
    res = 4 * (π / tan(s * (π / βsq)))
    res +=
        4 * sum(digamma_reg(-r + j / βsq) + digamma_reg(r + j / βsq) for j = (1-s):s)
    res -= sum(
        digamma_reg(
            1 // 2 +
            (lr == :left ? -1 : 1) *
            (Prs(r, j, β) + pm1 * V[a].dims[lr].P + pm2 * V[b].dims[lr].P) / β,
        ) for pm1 in (-1, 1) for pm2 in (-1, 1) for j = (1-s):2:(s-1) for
        (a, b) in ((1, 2), (3, 4)) for lr in (:left, :right)
    )
    res / β
end

function ell(V::NTuple{1,Field}, r, s)
    β = V[1].c.β
    βsq = β^2
    res = 4 * (π / tan(s * (π / β^2)))

    res +=
        4 * sum(digamma_reg(-r + j / βsq) + digamma_reg(r + j / βsq) for j = (1-s):s)
    res -=
        2 * sum(
            digamma_reg(
                1 // 2 +
                pm1 / β^2 / 2 +
                (pm2 * V[1].dims[lr].P + 2 * (lr == :left ? -1 : 1) * Prs(r, j, β)) /
                β,
            ) for j = (1-s):2:(s-1) for pm1 in (-1, 1) for pm2 in (-1, 1) for
            lr in (:left, :right)
        )
    res / β
end

function Base.show(io::IO, b::NCBlock)
    print(io, "G^($(b.chan))($(b.chan_field))")
end

function Base.show(io::IO, ::MIME"text/plain", b::NCBlock)
    print(io, "Non chiral")
    if typeof(b) <: LogBlock
        print(io, " logarithmic ")
    else
        print(io, " factorised ")
    end
    println(io, "block for the correlation")
    println(io, b.corr)
    println(io, "channel: $(b.chan), $(b.chan_field)")
end

#====================================================================================
Linear Combinations of blocks
====================================================================================#
function LinearCombinationBlock(bs::Vector{<:Block{T}}, coeffs) where {T}
    GenericLCBlock{T}(bs[1].chan_field, bs, coeffs, bs[1].chan)
end

function Base.:+(b1::LCBlock, b2::LCBlock)
    LCBlock(vcat(b1.blocks, b2.blocks), vcat(b1.coeffs, b2.coeffs))
end

function Base.:+(b1::LCBlock, b2::Block)
    LCBlock(vcat(b1.blocks, b2), vcat(b1.coeffs, 1))
end

function Base.:+(b1::Block, b2::LCBlock)
    LCBlock(vcat(b1, b2.blocks), vcat(1, b2.coeffs))
end

function Base.:+(b1::Block, b2::Block)
    LCBlock(vcat(b1, b2), vcat(1, 1))
end

function Base.:-(b1::Block, b2::Block)
    LCBlock(vcat(b1, b2), vcat(1, -1))
end

function Base.:*(a::Number, b::Block)
    LCBlock([b], [a])
end

function Base.show(io::IO, b::LCBlock)
    coeffs = [
        if c ≈ round(Int, real(c))
            round(Int, c)
        else
            c
        end for c in b.coeffs
    ]
    print(io, "($(coeffs[1])) * $(b.blocks[1])")
    for (i, bl) in enumerate(b.blocks[2:end])
        print(io, " + ($(coeffs[i+1])) * $bl")
    end
end
