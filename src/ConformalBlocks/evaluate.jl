struct PositionCache{T}
    x::T
    prefactor::T
    q::T
    logq::T
    q_powers::Vector{T}
end

struct LRPositionCache{T}
    left::PositionCache{T}
    right::PositionCache{T}
end

function PositionCache(x, ds::FourDimensions{T}, chan::Symbol, Nmax) where {T}
    ds = permute_dimensions(ds, chan)
    c = ds[1].c.c

    q = qfromx(x)
    
    a = (c-1)/24
    e0 = -ds[1].δ - ds[2].δ - a
    e1 = -ds[1].δ - ds[4].δ - a
    e2 = sum(ds[i].δ for i in 1:4) + a

    prefactor = x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2)
    if chan === :u
        prefactor *= x^(2ds[1].Δ)
    end

    sq = 16q
    q_powers = ones(T, Nmax+1)
    for i in 2:Nmax+1
        q_powers[i] = q_powers[i-1] * sq
    end

    return PositionCache{T}(x, prefactor, q, log(sq), q_powers)
end

function PositionCache(τ, ds::OneDimension{T}, chan::Symbol, Nmax) where {T}
    q = qfromτ(τ)
    prefactor = 1 / etaDedekind(complex(τ))
    q_powers = ones(T, Nmax+1)
    for i in 2:Nmax+1
        q_powers[i] = q_powers[i-1] * q
    end

    return PositionCache{T}(τ, prefactor, q, log(q), q_powers)
end

PositionCache(x, co::CorrelationChiral, chan) = PositionCache(x, co.dims, chan, co.Nmax)
PositionCache(x, co::CorrelationNonChiral, chan) = PositionCache(x, co[:left].dims, chan, co.Nmax)
PositionCache(x, b::Block) = PositionCache(x, b.corr, b.channel)

function LRPositionCache(x, co::Correlation{T}, chan) where {T}
    xbar = conj_q(x, co)
    return LRPositionCache{T}(
        PositionCache(x, co, chan), PositionCache(xbar, co, chan)
    )
end
LRPositionCache(x, b::Block) = LRPositionCache(x, b.corr, b.channel)

@inline function evalpoly(x::PositionCache, coeffs::Vector{T}) where {T}
    acc = zero(T)
    for i in 1:length(coeffs)
        acc += coeffs[i] * x.q_powers[i]
    end
    return acc
end

@inline evaluate_series(b::BlockChiral, x::PositionCache) = evalpoly(x, b._coefficients)
@inline evaluate_series_der(b::BlockChiral, x::PositionCache) = evalpoly(x, b._coefficients_der)
evaluate_series(b::BlockChiral, x::Number) = evaluate_series(b, PositionCache(x, b.corr, b.channel))
evaluate_series_der(b::BlockChiral, x::Number) = evaluate_series_der(b, PositionCache(x, b.corr, b.channel))

function crossratio(chan, x)
    chan === :s && return x
    chan === :t && return 1 - x
    chan === :u && return 1 / x
    error(
        """Incorrect channel specification in crossratio(channel, x):
        must be in $channels"""
    )
end
function modular_param(chan, τ)
    chan === :s && return τ
    chan === :t && return -1/τ
    chan === :u && return τ+1
end
channel_position(x, V::FourDimensions, chan) = crossratio(chan, x)
channel_position(x, V::OneDimension, chan) = modular_param(chan, x)
channel_position(x, c::CorrelationChiral, chan) = channel_position(x, c.dims, chan)
channel_position(x, c::CorrelationNonChiral, chan) = channel_position(x, c[:left].dims, chan)

conj_q(x, V::FourDimensions) = conj(x)
conj_q(τ, V::OneDimension) = -conj(τ)
conj_q(x, co::CorrelationChiral) = conj_q(x, co.dims)
conj_q(x, co::CorrelationNonChiral) = conj_q(x, co[:left].dims)

total_prefactor(b::BlockChiral, d::FourDimensions, x::PositionCache) = x.prefactor * (x.q_powers[2])^b.channel_dimension.δ
total_prefactor(b::BlockChiral, d::OneDimension, x::PositionCache) = x.prefactor * x.q^b.channel_dimension.δ
total_prefactor(b, x::PositionCache) = total_prefactor(b, b.corr.dims, x)
total_prefactor(b, x::Number) = total_prefactor(b, b.corr.dims, PositionCache(x, b))

function evaluate(b::BlockChiral{T}, x::PositionCache)::T where {T}
    d = b.channel_dimension
    p = total_prefactor(b, x)
    h = evaluate_series(b, x)

    if isdegenerate(d)
        # add log(q or 16q) * \sum C^N_rs (q or 16q)^N
        h += x.logq * evalpoly(x, b._missing_terms)
    end

    return p * h
end

function evaluate_der(b::BlockChiral{T}, x::PositionCache)::T where {T}
    d = b.channel_dimension
    qor16q = x.q_powers[2]
    p = x.prefactor * (qor16q)^b.channel_dimension.δ
    h = evaluate_series(b, x)
    hprime = evaluate_series_der(b, x)
    h = muladd(h, 2 * d.P * x.logq, hprime) # H_der = 2*P*log(q or 16q)*H + H'
    return p * h
end

@inline evaluate_lr(bs::LeftRight{BlockChiral}, x::LRPositionCache) =
    evaluate(bs[:left], x.left), evaluate(bs[:right], x.right)
@inline evaluate_lr_der(bs::LeftRight{BlockChiral}, x::LRPositionCache) =
    evaluate_der(bs[:left], x.left), evaluate_der(bs[:right], x.right)
@inline evaluate_lr(b::BlockFactorized, x) = evaluate_lr(b.chiral_blocks, x)
@inline evaluate_lr(b::BlockLogarithmic, x) = evaluate_lr(b.chiral_blocks, x)
@inline evaluate_lr_op(b::BlockLogarithmic, x) = evaluate_lr(b.chiral_blocks_op, x)
@inline evaluate_lr_der(b::BlockLogarithmic, x) = evaluate_lr_der(b.chiral_blocks_der, x)

@inline function evaluate(b::BlockFactorized{T}, x::LRPositionCache)::T where {T}
    return prod(evaluate_lr(b, x))
end

function evaluate(b::BlockLogarithmic{T}, x::LRPositionCache)::T where {T}
    V = b.channel_field
    r, s = V.indices
    s < 0 && return zero(T) # by convention G_(r, s<0) = 0
    Prs = V.P[:left]

    Freg, Fbar = evaluate_lr(b, x)
    F, Fregbar = evaluate_lr_op(b, x)

    if isaccidentallynonlogarithmic(b)
        Rreg = b.corr._Rmn_reg[:left][b.channel][r, s]
        Rregbar = b.corr._Rmn_reg[:right][b.channel][r, s]
        nbzeros = Rmn_zero_order(r, s, b.chiral_blocks[:left].dims)
        Freg * Fbar + (-1)^nbzeros * Rreg / Rregbar * F * Fregbar

    elseif islogarithmic(b)
        Fder, Fderbar = evaluate_lr_der(b, x)

        R = b.corr._Rmn[:left][b.channel][r, s]
        Rbar = b.corr._Rmn[:right][b.channel][r, s]
        l = b.ell

        (Freg - R / 2 / Prs * Fder) * Fbar +
        R / Rbar * F * (Fregbar - Rbar / 2 / Prs * Fderbar) +
        R / 2 / Prs * l * F * Fbar
    end
end

function evaluate(b::BlockInterchiral{T}, x)::T where {T}
    res = zero(T)
    for i in eachindex(b.blocks)
        res += evaluate((b.blocks)[i], x) .* b.shifts[i]
    end
    return res
end

evaluate(b::BlockChiral, x::Number) = evaluate(b, PositionCache(x, b))
evaluate(b::BlockNonChiral, x::Number) = evaluate(b, LRPositionCache(x, b))
evaluate_der(b::BlockChiral, x::Number) = evaluate_der(b, PositionCache(x, b))
evaluate_der(b::BlockNonChiral, x::Number) = evaluate_der(b, LRPositionCache(x, b))
