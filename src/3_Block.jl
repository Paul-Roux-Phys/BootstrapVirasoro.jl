"""
Abstract supertype for all conformal block types.

# Hierarchy

```text
Block  (Abstract Type)
├─ ChiralBlock
├─ NonChiralBlock  (Abstract Type)
│  ├─ Factorised block
│  └─ LogBlock
└─ LinearCombinationBlock
```
"""
abstract type Block end # general conformal block. Can be interchiral, non-chiral or chiral

"""
# Type

Type to represent the series expansion of a chiral block.
Aliased to CBlock.

# Constructors 

        ChiralBlock(::ChiralCorrelation4, ::Symbol, ::CD, Δmax=co.Δmax::Int, der=false) # 4pt 𝕊² chiral block
        ChiralBlock(::ChiralCorrelation1, ::CD, Δmax=co.Δmax::Int, der=false)           # 1pt 𝕋² chiral block

Compute the series coefficients of a chiral block associated to the
 [`ChiralCorrelation`](@ref Correlation) `co`, in the channel `chan`.
If `V` is a degenerate field, compute the ``P``-regularisation of the block instead.
If `der=true`, also compute the coefficients of the series expansion of the derivative of the block.

Aliased to CBlock.

## Arguments

- A chiral correlation object
- A channel for the block. Only required for four-point blocks.
- A channel `ConformalDimension` 
- `Δmax`: integer up to which the series is evaluated.
- `der`: whether to compute the coefficients of the block's derivative as well.
"""
struct ChiralBlock <: Block
    corr::CCo
    chan_field::CD
    coeffs::Vector{Acb}
    coeffs_der::Vector{Acb}
    missing_terms::Vector{Acb}
    chan::Symbol
    Δmax::Int
end
const CBlock = ChiralBlock

"""
Holds all of the precomputable data about a position that is needed
to evaluate blocks:
powers and log of the nome ``q`` or of ``16q``, value of the part of
the prefactor of the block that is independent
of the channel dimension.
"""
struct ChiralPosCache
    x::Acb
    prefactor::Acb
    q::Acb
    logq::Acb
    q_powers::Vector{Acb}
end

"""
# Type

Abstract supertype to represent the series expansion of a non chiral block,
either left-right factorised or logarithmic.

Aliased to NCBlock.

# Constructors

        NCBlock(::NonChiralCorrelation, ::Symbol, ::Field, Δmax=co.Δmax::Int) # 4pt 𝕊² chiral block
        NCBlock(::NonChiralCorrelation, ::Field, Δmax=co.Δmax::Int)           # 1pt 𝕋² chiral block

Compute the series coefficients of the left and right blocks associated to the
 non chiral correlation `co`, in the channel `chan`.
If the channel field is degenerate, compute the ``P``-regularisation of the block instead.
If the channel field is a non-diagonal field of the type ``V_{(r, s>0)}``,
the block is logarithmic, except if the residue ``R_{r,s}`` vanishes.

## Arguments

- A `NonChiralCorrelation` object
- A channel for the block. Only required for four-point blocks.
- A channel `Field` 
- `Δmax`: integer up to which the series is evaluated. Defaults to the correlation's `Δmax`.
"""
abstract type NonChiralBlock <: Block end
const NCBlock = NonChiralBlock

struct FactorizedBlock <: NonChiralBlock
    chan_field::Field
    corr::Corr
    cblocks::LR{CBlock}
    Δmax::Int
    chan::Symbol
end

struct LogBlock <: NonChiralBlock
    cblocks::LR{CBlock} # blocks for V_(r, s) x V_(r, -s)
    cblocks_op::LR{CBlock} # blocks for V_(r, -s) x V_(r, s)
    cblocks_der::Union{LR{CBlock},Nothing}
    corr::Corr
    chan_field::Field
    R::Acb
    Rbar::Acb
    ell::Acb
    Δmax::Int
    nbzeros::Int
    chan::Symbol
end

const NonChiralPosCache = LeftRight{ChiralPosCache}

"""
# Type

Type to represent linear combinations of blocks. Contains

- an array of `Block`s
- an array of coefficients of the linear combination.

# Constructors

`LinearCombinationBlock`s are constructed by summing blocks

## Examples

```julia
# define a correlation co, two fields V1, V2.
b1 = CBlock(co, :s, V1)
b2 = CBlock(co, :s, V2)
b = b1 - 2b2 # LinearCombinationBlock
```
"""
abstract type LinearCombinationBlock <: Block end
const LCBlock = LinearCombinationBlock

struct GenericLCBlock <: LCBlock
    chan_field::Field
    blocks::Vector{Block}
    coeffs::Vector{Acb}
    chan::Symbol
end

qfromx(x) = exp(-(π * ellipticK(1 - x) / ellipticK(x)))
xfromq(q) = jtheta2(0, q)^4 / jtheta3(0, q)^4
qfromτ(τ) = exp(2 * im * (π * τ))
τfromx(x) = (log(qfromx(x)) / π) / im
xfromτ(τ) = xfromq(exp(im * (π * τ)))

"""
        (b::Block)(x)
        
Evaluate the block `b` at the parameter `x`.
If the block is a four-point block, `x` is the cross-ratio.
If the block is a one-point block, `x` is the torus' modulus.
"""
function (b::Block)(x) end

#=============================================================================
Chiral Blocks
=============================================================================#
function series_H(d::CD, Δmax, CNmn) 
    P = d.P
    P2 = P^2
    isKac = d.isKac
    r, s = d.r, d.s

    coeffs = [zero(Acb) for _ = 1:Δmax+2, _ = 1:Δmax+2]
    all_mns = union([CNmn.keys[N] for N = 1:Δmax]...)
    buf = zero(Acb)
    o = one(Acb)
    four = convert(Acb, 4)
    for (m, n) in all_mns
        if isKac && m == r && n == s
            # coeffs[m, n] = -inv(4CNmn.δs[r, s])
            buf = MA.operate_to!!(buf, *, four, CNmn.δs[r, s]) # 4 * δrs
            buf = MA.operate_to!!(buf, /, o, buf) # 1/ (4*δrs)
            buf = MA.operate!!(-, buf) # -1/ (4*δrs)
            MA.operate_to!!(coeffs[m, n], copy, buf) # store the result
        else
            # coeffs[m, n] = 1 / (P2 - CNmn.δs[m, n])
            buf = MA.operate_to!!(buf, -, P2, CNmn.δs[m, n])
            buf = MA.operate_to!!(buf, /, o, buf)
            coeffs[m, n] = MA.operate_to!!(coeffs[m, n], copy, buf)
        end
    end

    H = [zero(Acb) for _ = 1:Δmax+1]
    H[1] = one(Acb)
    buf = zero(Acb)
    for N = 1:Δmax
        for (m, n) in CNmn.keys[N]
            (N > Δmax || m > Δmax || n > Δmax) && continue
            # H[N+1] += CNmn[N, m, n] * coeffs[m, n]
            buf = MA.operate_to!!(buf, *, CNmn[N, m, n], coeffs[m, n])
            buf = MA.operate_to!!(buf, +, H[N+1], buf)
            H[N+1] = MA.operate_to!!(H[N+1], copy, buf)
        end
    end

    return H
end

function series_H_der(d::CD, Δmax, CNmn) 
    P = d.P
    P2 = P^2
    mtwoP = -2P
    buf = zero(Acb)

    H = [zero(Acb) for _ = 1:Δmax+1]

    # coeffs = Dict(
    #     (m, n) => -2P / (P2 - CNmn.δs[m, n])^2 for
    #     (m, n) in union([CNmn.keys[N] for N = 1:Δmax]...)
    # )

    all_mns = union([CNmn.keys[N] for N = 1:Δmax]...)
    coeffs = Dict((m, n) => zero(Acb) for (m, n) in all_mns)
    for (m, n) in all_mns
        # -2P / (P2 - CNmn.deltas[m, n])^2
        buf = MA.operate_to!!(buf, -, P2, CNmn.δs[m, n])
        buf = MA.operate_to!!(buf, *, buf, buf)
        buf = MA.operate_to!!(buf, /, mtwoP, buf)
        coeffs[(m, n)] = MA.operate_to!!(coeffs[(m, n)], copy, buf)
    end

    for N = 1:Δmax
        for (m, n) in CNmn.keys[N]
            (N > Δmax || m > Δmax || n > Δmax) && continue
            # H[N+1] += CNmn[N, m, n] * coeffs[(m, n)]
            buf = MA.operate_to!!(buf, *, CNmn[N, m, n], coeffs[(m, n)])
            buf = MA.operate_to!!(buf, +, H[N+1], buf)
            H[N+1] = MA.operate_to!!(H[N+1], copy, buf)
        end
    end

    return H
end

function ChiralBlock(
    co::CCo,
    chan::Symbol,
    d::CD,
    Δmax = missing,
    der = false,
) 
    Δmax === missing && (Δmax = co.Δmax)
    CNmn = getCNmn(co, chan)
    coeffs = series_H(d, Δmax, CNmn)
    coeffs_der = Vector{Acb}()
    der && (coeffs_der = series_H_der(d, Δmax, CNmn))
    if !d.degenerate
        missing_terms = Vector{Acb}()
    else
        r, s = d.r, d.s
        missing_terms = [
            (N > 0 && (r, s) in CNmn.keys[N]) ? CNmn[N, r, s] : zero(Acb) for N = 0:Δmax
        ]
    end

    CBlock(co, d, coeffs, coeffs_der, missing_terms, chan, Δmax)
end

function ChiralBlock(
    co::ChiralCorrelation1,
    d::CD = CD(),
    Δmax = missing,
    der = false,
)
    CBlock(co, :τ, d, Δmax, der)
end

getRmn(b::CBlock) = getRmn(b.corr, b.chan)
getRmnreg(b::CBlock) = getRmnreg(b.corr, b.chan)
getCNmn(b::CBlock) = getCNmn(b.corr, b.chan)
getc(b::CBlock) = b.corr.c

function ChiralPosCache(x, ds::NTuple{4,CD}, chan::Symbol, Δmax) 
    ds = permute_4(ds, chan)

    q = qfromx(x)

    e0 = -ds[1].Δ - ds[2].δ
    chan === :u && (e0 += 2ds[1].Δ)
    e1 = -ds[1].Δ - ds[4].δ
    e2 = sum(ds[i].δ for i = 1:3) + ds[4].Δ

    prefactor = x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2)

    sq = 16q
    q_powers = ones(Acb, Δmax + 1)
    for i = 2:(Δmax+1)
        q_powers[i] = q_powers[i-1] * sq
    end

    return ChiralPosCache(x, prefactor, q, log(sq), q_powers)
end

function ChiralPosCache(τ, _::NTuple{1,CD}, _::Symbol, Δmax) 
    q = qfromτ(τ)
    prefactor = 1 / etaDedekind(complex(τ))
    q_powers = ones(Acb, Δmax + 1)
    for i = 2:(Δmax+1)
        q_powers[i] = q_powers[i-1] * q
    end

    return ChiralPosCache(τ, prefactor, q, log(q), q_powers)
end

ChiralPosCache(x, co::CCo, chan) = ChiralPosCache(x, co.fields, chan, co.Δmax)
PosCache(x, co::CCo, chan) = ChiralPosCache(x, co.fields, chan, co.Δmax)
PosCache(x, b::CBlock) = PosCache(x, b.corr, b.chan)

function evalpoly(x::ChiralPosCache, coeffs::Vector{Acb}) 
    res = zero(Acb)
    buf = zero(Acb)
    for i = 1:length(coeffs)
        # res += coeffs[i] * x.q_powers[i]
        buf = MA.operate_to!!(buf, *, coeffs[i], x.q_powers[i])
        res = MA.operate!!(+, res, buf)
    end
    return res
end

eval_series(b::CBlock, x::ChiralPosCache) = evalpoly(x, b.coeffs)
eval_series_der(b::CBlock, x::ChiralPosCache) = evalpoly(x, b.coeffs_der)
eval_series(b::CBlock, x::Number) = eval_series(b, PosCache(x, b.corr, b.chan))
eval_series_der(b::CBlock, x::Number) =
    eval_series_der(b, PosCache(x, b.corr, b.chan))

prefactor(b::CBlock, x::Number) = PosCache(x, b).prefactor
total_prefactor(b::CBlock, x::ChiralPosCache, _::Correlation4) =
    x.prefactor * (x.q_powers[2])^b.chan_field.δ
total_prefactor(b::CBlock, x::ChiralPosCache, _::Correlation1) =
    x.prefactor * x.q^b.chan_field.δ
total_prefactor(b::CBlock, x::ChiralPosCache) = total_prefactor(b, x, b.corr)
total_prefactor(b::CBlock, x::Number) = total_prefactor(b, PosCache(x, b))

function (b::CBlock)(x::ChiralPosCache)::Acb 
    d = b.chan_field
    p = total_prefactor(b, x)
    h = eval_series(b, x)

    if d.degenerate
        # add log(q or 16q) * \sum C^N_rs (q or 16q)^N
        h += x.logq * evalpoly(x, b.missing_terms)
    end

    return p * h
end

function (b::CBlock)(x::ChiralPosCache, _::Bool)::Acb 
    d = b.chan_field
    qor16q = x.q_powers[2]
    p = x.prefactor * (qor16q)^b.chan_field.δ
    h = eval_series(b, x)
    hprime = eval_series_der(b, x)
    h = muladd(h, 2 * d.P * x.logq, hprime) # H_der = 2*P*log(q or 16q)*H + H'
    return p * h
end

(b::CBlock)(x::Number) = b(PosCache(x, b))
(b::CBlock)(x::Number, _::Bool) = b(PosCache(x, b), true)

function Base.show(io::IO, ::MIME"text/plain", b::CBlock)
    println(io, "Chiral block for the correlation")
    println(io, b.corr)
    println(io, "Channel: $(b.chan), $(b.chan_field)")
end

function Base.show(io::IO, b::CBlock)
    println(io, "F^($(b.chan))($(b.chan_field))")
end

#================================================================================
Non Chiral blocks
================================================================================#
function FactorizedBlock(co::NCCo, chan, V, Δmax) 
    bl = CBlock(co[:left], chan, V[:left], Δmax)
    br = CBlock(co[:right], chan, V[:right], Δmax)
    FactorizedBlock(V, co, LR(bl, br), Δmax, chan)
end

# function islogarithmic(V::Field)
# r, s = V.r, V.s
# !V.diagonal && r * s != 0 && r % 1 == s % 1 == 0 && r > 0 && s > 0
# end

islogarithmic(V::Field) = !V.diagonal && V[:left].degenerate
islogarithmic(b::NCBlock) = islogarithmic(b.chan_field)

function isaccidentallynonlogarithmic(co::NCCo, chan, V)
    !islogarithmic(V) && return false
    return Rmn_zero_order(V.r, V.s, getfields(co[:left], chan)) > 0
end

isaccidentallynonlogarithmic(b::NCBlock) =
    isaccidentallynonlogarithmic(b.corr, b.chan, b.chan_field)

function LogBlock(co::NCCo, chan, V, Δmax) 
    V_op = swap_lr(V)
    VV = V.s > 0 ? (V, V_op) : (V_op, V) # (V_(r, s > 0), V_(r, -s))
    r, s = VV[1].r, VV[1].s
    left, left_op, right, right_op = Tuple(
        CBlock(co[lr], chan, v.dims[lr], Δmax, false) for lr in (:left, :right)
        for v in VV
    )
    R, Rbar, l = zero(Acb), zero(Acb), zero(Acb)
    cl, cr = co[:left], co[:right]
    if isaccidentallynonlogarithmic(co, chan, V)
        chiral_blocks_der = nothing
        if V.s > 0
            R = getRmnreg(cl, chan)[r, s]
            Rbar = getRmnreg(cr, chan)[r, s]
        end
        nbzeros = Rmn_zero_order(r, s, left.corr.fields)
        l = zero(Acb)
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
    LogBlock(
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

function NonChiralBlock() end

function NonChiralBlock(co::NCCo, chan, V, Δmax = missing)
    Δmax === missing && (Δmax = co.Δmax)
    if islogarithmic(V)
        LogBlock(co, chan, V, Δmax)
    else
        FactorizedBlock(co, chan, V, Δmax)
    end
end
NonChiralBlock(co::NonChiralCorrelation1, V, Δmax = missing) =
    NCBlock(co, :τ, V, Δmax)

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

function LeftRight{ChiralPosCache}(x, co::NCCo, chan) 
    xbar = conj_q(x, co)
    return NonChiralPosCache(PosCache(x, co[:left], chan), PosCache(xbar, co[:right], chan))
end
LeftRight{ChiralPosCache}(x, b::NCBlock) = NonChiralPosCache(x, b.corr, b.chan)

PosCache(x, co::NCCo, chan) = NonChiralPosCache(x, co, chan)
PosCache(x, b::NCBlock) = NonChiralPosCache(x, b)

conj_q(x, _::Correlation4) = conj(x)
conj_q(τ, _::Correlation1) = -conj(τ)

prefactor(b::NonChiralBlock, x::NonChiralPosCache) =
    prefactor(b.cblocks.left, x.left) * prefactor(b.cblocks.right, x.right)
prefactor(b::NonChiralBlock, x::Number) = prefactor(b, NonChiralPosCache(x, b))

eval_lr(bs::LR{CBlock}, x)  = bs.left(x.left), bs.right(x.right)
eval_lr_der(bs::LeftRight{CBlock}, x)  =
    bs.left(x.left, true), bs.right(x.right, true)
eval_lr(b::FactorizedBlock, x)  = eval_lr(b.cblocks, x)
eval_lr(b::LogBlock, x) = eval_lr(b.cblocks, x)
eval_lr_op(b::LogBlock, x) = eval_lr(b.cblocks_op, x)
eval_lr_der(b::LogBlock, x) = eval_lr_der(b.cblocks_der, x)

function (b::FactorizedBlock)(x::NonChiralPosCache)::Acb 
    lr = eval_lr(b, x)
    return lr[1] * lr[2]
end

function (b::LogBlock)(x::NonChiralPosCache)::Acb 
    V = b.chan_field
    s = V.s
    s < 0 && return zero(Acb) # by convention G_(r, s<0) = 0
    Prs = V[:left].P

    Freg, Fbar = eval_lr(b, x)
    F, Fregbar = eval_lr_op(b, x)
    R, Rbar = b.R, b.Rbar

    if isaccidentallynonlogarithmic(b)
        Freg * Fbar + R / Rbar * F * Fregbar

    elseif islogarithmic(b)
        Fder, Fderbar = eval_lr_der(b, x)

        l = b.ell

        (Freg - R / 2 / Prs * Fder) * Fbar +
        R / Rbar * F * (Fregbar - Rbar / 2 / Prs * Fderbar) +
        R / 2 / Prs * l * F * Fbar
    end
end

(b::NCBlock)(x::Number) = b(NonChiralPosCache(x, b))
(b::NCBlock)(x::Number, _::Bool) = b(NonChiralPosCache(x, b), true)

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
function LinearCombinationBlock(bs::Vector{<:Block}, coeffs) 
    GenericLCBlock(bs[1].chan_field, bs, coeffs, bs[1].chan)
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

function (b::LCBlock)(x)::Acb 
    res = zero(Acb)
    for i in eachindex(b.blocks)
        res += (b.blocks)[i](x) .* b.coeffs[i]
    end
    return res
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
