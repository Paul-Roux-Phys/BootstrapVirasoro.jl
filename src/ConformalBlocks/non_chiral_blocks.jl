abstract type NonChiralBlock{T,U} <: Block{T,U} end

struct FactorizedBlock{T,U} <: NonChiralBlock{T,U}

    Nmax::Int
    correlation::Correlation
    channel_field::Field{T}
    chiral_blocks::LeftRight{ChiralBlock{T,U}}

end

struct LogarithmicBlock{T,U} <: NonChiralBlock{T,U}

    Nmax::Int
    correlation::Correlation
    channel_field::Field{T}
    chiral_blocks::LeftRight{ChiralBlock{T,U}} # blocks for V_(r, s) x V_(r, -s)
    chiral_blocks_op::LeftRight{ChiralBlock{T,U}} # blocks for V_(r, -s) x V_(r, s)
    chiral_blocks_der::Union{LeftRight{ChiralBlock{T,U}},Nothing}
    R::T
    Rbar::T
    ell::T
    nbzeros::Int

end

function FactorizedBlock(co::Correlation{T,U}, chan, V, Nmax) where {T,U}
    left, right = Tuple(ChiralBlock(co, chan, V, lr, Nmax) for lr in (:left, :right))
    W = typeof(co[:left]).parameters[2]
    FactorizedBlock{T,W}(Nmax, co, V, LeftRight((left, right)))
end

function islogarithmic(V::Field)
    r, s = indices(V)
    !isdiagonal(V) && r*s != 0 && r%1 == s%1 == 0 && r > 0 && s > 0
end

islogarithmic(b::Block) = islogarithmic(b.channel_field)

function isaccidentallynonlogarithmic(co, chan, V)
    !islogarithmic(V) && return false
    return Rmn_zero_order(V.r, V.s, permute_dimensions(co[:left].fields, chan)) > 0
end

isaccidentallynonlogarithmic(b) =
    isaccidentallynonlogarithmic(b.corr, b.channel, b.channel_field)

function LogarithmicBlock(co::CorrelationNonChiral{T,U}, chan, V, Nmax) where {T,U}
    V_op = swap_lr(V)
    VV = V.s > 0 ? (V, V_op) : (V_op, V) # (V_(r, s > 0), V_(r, -s))
    r, s = VV[1].r, VV[1].s
    left, left_op, right, right_op =
        Tuple(ChiralBlock(co, chan, v, lr, Nmax) for lr in (:left, :right) for v in VV)
    R, Rbar, l = zero(T), zero(T), zero(T)
    if isaccidentallynonlogarithmic(co, chan, V)
        chiral_blocks_der = nothing
        if V.s > 0
            R = co._Rmn_reg[:left][chan][r, s]
            Rbar = co._Rmn_reg[:right][chan][r, s]
        end
        nbzeros = Rmn_zero_order(r, s, left.fields)
        l = 0
    else
        leftder = ChiralBlock(co, chan, VV[2], :left, Nmax, der = true)
        rightder = ChiralBlock(co, chan, VV[1], :right, Nmax, der = true)
        chiral_blocks_der = LeftRight((leftder, rightder))
        if V.s > 0
            R = co._Rmn[:left][chan][r, s]
            Rbar = co._Rmn[:right][chan][r, s]
            l = ell(co, chan, r, VV[1].s)
        end
        nbzeros = 0
    end
    W = typeof(co[:left]).parameters[2]
    LogarithmicBlock{T,W}(
        Nmax,
        co,
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

reflect(b::NonChiralBlock) = FactorizedBlock(b.corr, b.channel, reflect(b.channel_field), b.Nmax)

function NonChiralBlock(co, chan, V, Nmax)
    if islogarithmic(V)
        LogarithmicBlock(co, chan, V, Nmax)
    else
        FactorizedBlock(co, chan, V, Nmax)
    end
end

function Base.getproperty(b::NonChiralBlock, s::Symbol)
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

function Base.getindex(b::NonChiralBlock, s::Symbol)
    s in (:left, :right) && return getfield(b, :chiral_blocks)[s]
    return error("invalid index")
end

function Base.show(io::IO, ::MIME"text/plain", b::NonChiralBlock)
    print(io, "Non chiral")
    if typeof(b) <: LogarithmicBlock
        print(io, " logarithmic ")
    else
        print(io, " factorised ")
    end
    println(io, "block for the correlation")
    println(io, b.corr)
    println(io, "channel: $(b.channel), $(b.channel_field)")
end

function Base.show(io::IO, b::NonChiralBlock)
    print(io, "G^($(b.channel))($(b.channel_field))")
end

"""Factor ``\\ell_{(r,s)}`` that appears in logarithmic blocks"""
ell(b::Block, r, s) = ell(b.corr.fields, b.channel, r, s)

function ell(V::FourFields, chan, r, s)
    V = permute_fields(V, chan)
    β = V[1].c.β
    βsq = β^2
    res = 4 * (π / tan(s * (π / β^2)))
    res +=
        4 * sum(
        digamma_reg(-r + j/βsq) + digamma_reg(r + j/βsq) for
        j = (1-s):s
    )
    res -= sum(
        digamma_reg(
            1 // 2 + (lr == :left ? -1 : 1) *
                (P_rs(r, j, β) + pm1 * V[a].dims[lr].P + pm2 * V[b].dims[lr].P) / β
        )
        for pm1 in (-1, 1) for pm2 in (-1, 1) for j = (1-s):2:(s-1) for
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
                (pm2 * V[1].dims[lr].P + 2 * (lr == :left ? -1 : 1) * P_rs(r, j, β)) / β,
            ) for j = (1-s):2:(s-1) for pm1 in (-1, 1) for pm2 in (-1, 1) for
            lr in (:left, :right)
        )
    res / β
end

ell(c::Correlation, chan, r, s) = ell(c.fields, chan, r, s)
ell(b::Block) = ell(b.corr, b.channel, b.channel_field.r, b.channel_field.s)
