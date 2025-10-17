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

function FactorizedBlock(co::NCCo{T}, chan, V, Δmax) where {T}
    bl = Block(co[:left], chan, V.dims.left, Δmax)
    br = Block(co[:right], chan, V.dims.right, Δmax)
    FactorizedBlock{T}(V, co, LR(bl, br), Δmax, chan)
end

function islogarithmic(V::Field)
    r, s = indices(V)
    !V.diagonal && r * s != 0 && r % 1 == s % 1 == 0 && r > 0 && s > 0
end

islogarithmic(b::NCBlock) = islogarithmic(b.chan_field)

function isaccidentallynonlogarithmic(co::NCCo, chan, V)
    !islogarithmic(V) && return false
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
        nbzeros = Rmn_zero_order(r, s, left.fields)
        l = zero(T)
    else
        leftder = CBlock(cl, chan, VV[2].dims.left, Δmax, true)
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

function NonChiralBlock(co, chan, V, Δmax)
    if islogarithmic(V)
        LogBlock(co, chan, V, Δmax)
    else
        FactorizedBlock(co, chan, V, Δmax)
    end
end

Base.getindex(b::NCBlock, s::Symbol) = b.cblocks[s]

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

function Base.show(io::IO, b::NCBlock)
    print(io, "G^($(b.chan))($(b.chan_field))")
end

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
        4 * sum(
    digamma_reg(-r + j/βsq) + digamma_reg(r + j/βsq)
            for j = (1-s):s
        )
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
