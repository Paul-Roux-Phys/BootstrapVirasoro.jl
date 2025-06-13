abstract type BlockNonChiral{T,U} <: Block{T,U} end

struct BlockFactorized{T,U} <: BlockNonChiral{T,U}

    Nmax::Int
    channel_field::Field{T}
    chiral_blocks::LeftRight{BlockChiral{T,U}}

end

struct BlockLogarithmic{T,U} <: BlockNonChiral{T,U}

    Nmax::Int
    channel_field::Field{T}
    chiral_blocks::LeftRight{BlockChiral{T,U}} # blocks for V_(r, s) x V_(r, -s)
    chiral_blocks_op::LeftRight{BlockChiral{T,U}} # blocks for V_(r, -s) x V_(r, s)
    chiral_blocks_der::Union{LeftRight{BlockChiral{T,U}},Nothing}
    R::T
    Rbar::T
    ell::T
    nbzeros::Int

end

function BlockFactorized(co::Correlation{T,U}, chan, V, Nmax) where {T,U}
    left, right = Tuple(BlockChiral(co, chan, V, lr, Nmax) for lr in (:left, :right))
    W = typeof(co[:left]).parameters[2]
    BlockFactorized{T,W}(Nmax, V, LeftRight((left, right)))
end

function islogarithmic(V::Field)
    r, s = V.indices
    !isdiagonal(V) && r*s != 0 && r%1 == s%1 == 0
end

islogarithmic(b::Block) = islogarithmic(b.channel_field)

function isaccidentallynonlogarithmic(co, chan, V)
    !islogarithmic(V) && return false
    return Rmn_zero_order(V.r, V.s, permute_dimensions(co[:left].dims, chan)) > 0
end

isaccidentallynonlogarithmic(b) =
    isaccidentallynonlogarithmic(b.corr, b.channel, b.channel_field)

function BlockLogarithmic(co::CorrelationNonChiral{T,U}, chan, V, Nmax) where {T,U}
    V_op = swap_lr(V)
    VV = V.s > 0 ? (V, V_op) : (V_op, V) # (V_(r, s > 0), V_(r, -s))
    r, s = VV[1].r, VV[1].s
    left, left_op, right, right_op =
        Tuple(BlockChiral(co, chan, v, lr, Nmax) for lr in (:left, :right) for v in VV)
    R, Rbar, l = zero(T), zero(T), zero(T)
    if isaccidentallynonlogarithmic(co, chan, V)
        chiral_blocks_der = nothing
        if V.s > 0
            R = co._Rmn_reg[:left][chan][r, s]
            Rbar = co._Rmn_reg[:right][chan][r, s]
        end
        nbzeros = Rmn_zero_order(r, s, left.dims)
        l = 0
    else
        leftder = BlockChiral(co, chan, VV[2], :left, Nmax, der = true)
        rightder = BlockChiral(co, chan, VV[1], :right, Nmax, der = true)
        chiral_blocks_der = LeftRight((leftder, rightder))
        if V.s > 0
            R = co._Rmn[:left][chan][r, s]
            Rbar = co._Rmn[:right][chan][r, s]
            l = ell(co, chan, r, VV[1].s)
        end
        nbzeros = 0
    end
    W = typeof(co[:left]).parameters[2]
    BlockLogarithmic{T,W}(
        Nmax,
        V,
        LeftRight((left, right)),
        LeftRight((left_op, right_op)),
        chiral_blocks_der,
        R,
        Rbar,
        l,
        nbzeros,
    )
end

function BlockNonChiral(co, chan, V, Nmax)
    if islogarithmic(V)
        BlockLogarithmic(co, chan, V, Nmax)
    else
        BlockFactorized(co, chan, V, Nmax)
    end
end

function Base.getproperty(b::BlockNonChiral, s::Symbol)
    s === :c && return getfield(b.channel_field.dims[:left], s)
    s === :corr && begin
        left_corr = getfield(b, :chiral_blocks)[:left].corr
        right_corr = getfield(b, :chiral_blocks)[:right].corr
        return Correlation(left_corr, right_corr)
    end
    s === :fields && return getproperty(b, :corr).fields
    s === :channel && return getproperty(getfield(b, :chiral_blocks)[:left], s)
    s in (:r, :s) && return getproperty(getfield(b, :channel_field), s)
    return getfield(b, s)
end

function Base.getindex(b::BlockNonChiral, s::Symbol)
    s in (:left, :right) && return getfield(b, :chiral_blocks)[s]
    return error("invalid index")
end

function Base.show(io::IO, ::MIME"text/plain", b::BlockNonChiral)
    print(io, "Non chiral")
    if typeof(b) <: BlockLogarithmic
        print(io, " logarithmic ")
    else
        print(io, " factorised ")
    end
    println(io, "block for the correlation")
    println(io, b.corr)
    println(io, "channel: $(b.channel), $(b.channel_field)")
end

function Base.show(io::IO, b::BlockNonChiral)
    print(io, "G^($(b.channel))($(b.channel_field))")
end

"""Factor ``\\ell_{(r,s)}`` that appears in logarithmic blocks"""
ell(b::Block, r, s) = ell(b.corr.fields, b.channel, r, s)

"""Factor \ell_{(r,s)} that appears in logarithmic blocks"""
function ell(V::FourFields, chan, r, s)
    V = permute_fields(V, chan)
    β = V[1].c.β
    res = 4 * (π / tan(s * (π / β^2)))
    res +=
        4 * sum(
            digamma_reg(-2 * P_rs(r, j, β) / β) + digamma_reg(2 * P_rs(r, -j, β) / β) for
            j = (1-s):s
        )
    res -= sum(
        digamma_reg(
            1 // 2 +
            (lr == :left ? -1 : 1) * (P_rs(r, j, β) + pm1 * V[a].P[lr] + pm2 * V[b].P[lr]) / β,
        ) for pm1 in (-1, 1) for pm2 in (-1, 1) for j = (1-s):2:(s-1) for
        (a, b) in ((1, 2), (3, 4)) for lr in (:left, :right)
    )
    res / β
end

function ell(V::OneField, chan, r, s)
    β = V[1].c.β
    res = 4 * (π / tan(s * (π / β^2)))
    res +=
        4 * sum(
            digamma_reg(-2 * P_rs(r, j, β) / β) + digamma_reg(2 * P_rs(r, -j, β) / β) for
            j = (1-s):s
        )
    res -=
        2 * sum(
            digamma_reg(
                1//2 +
                pm1 * inv(β)^2 / 2 +
                (pm2 * V[1].P[lr] + 2 * (lr == :left ? -1 : 1) * P_rs(r, j, β)) / β,
            ) for j = (1-s):2:(s-1) for pm1 in (-1, 1) for pm2 in (-1, 1) for
            lr in (:left, :right)
        )
    res / β
end

ell(c::Correlation, chan, r, s) = ell(c.fields, chan, r, s)
ell(b::Block) = ell(b.corr, b.channel, b.channel_field.r, b.channel_field.s)
