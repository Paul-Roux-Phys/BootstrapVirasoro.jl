"""Nome `q` from the cross-ratio `x`"""
@inline qfromx(x) = exp(- (π * ellipticK(1 - x) / ellipticK(x)))
"""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0, q)^4 / jtheta3(0, q)^4
@inline qfromτ(τ) = exp(im*(π*τ))

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
    e2 = sum(ds[i].δ for i = 1:4) + a

    prefactor = x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2)
    if chan === :u
        prefactor *= x^(2ds[1].Δ)
    end

    sq = 16q
    q_powers = ones(T, Nmax+1)
    for i = 2:(Nmax+1)
        q_powers[i] = q_powers[i-1] * sq
    end

    return PositionCache{T}(x, prefactor, q, log(sq), q_powers)
end

function PositionCache(τ, _::OneDimension{T}, _::Symbol, Nmax) where {T}
    q = qfromτ(τ)
    prefactor = 1 / etaDedekind(complex(τ))
    q_powers = ones(T, Nmax+1)
    for i = 2:(Nmax+1)
        q_powers[i] = q_powers[i-1] * q
    end

    return PositionCache{T}(τ, prefactor, q, log(q), q_powers)
end

PositionCache(x, co::CorrelationChiral, chan) = PositionCache(x, co.dims, chan, co.Nmax)
PositionCache(x, co::CorrelationNonChiral, chan) =
    PositionCache(x, co[:left].dims, chan, co.Nmax)
PositionCache(x, b::Block) = PositionCache(x, b.corr, b.channel)

function LRPositionCache(x, co::Correlation{T,U}, chan) where {T,U}
    xbar = conj_q(x, co)
    return LRPositionCache{T}(PositionCache(x, co, chan), PositionCache(xbar, co, chan))
end
LRPositionCache(x, b::Block) = LRPositionCache(x, b.corr, b.channel)

@inline function evalpoly(x::PositionCache, coeffs::Vector{T}) where {T}
    res = zero(T)
    for i = 1:length(coeffs)
        res += coeffs[i] * x.q_powers[i]
    end
    return res
end

@inline evaluate_series(b::ChiralBlock, x::PositionCache) = evalpoly(x, b._coefficients)
@inline evaluate_series_der(b::ChiralBlock, x::PositionCache) =
    evalpoly(x, b._coefficients_der)
evaluate_series(b::ChiralBlock, x::Number) =
    evaluate_series(b, PositionCache(x, b.corr, b.channel))
evaluate_series_der(b::ChiralBlock, x::Number) =
    evaluate_series_der(b, PositionCache(x, b.corr, b.channel))

function crossratio(chan, x)
    chan === :s && return x
    chan === :t && return 1 - x
    chan === :u && return 1 / x
    error("""Incorrect channel specification in crossratio(channel, x):
          must be in $channels""")
end
function modular_param(chan, τ)
    chan === :s && return τ
    chan === :t && return -1/τ
    chan === :u && return τ+1
end
channel_position(x, _::Correlation{T,U}, chan,) where {T,U<:FourPoints} =
    crossratio(chan, x)
channel_position(x, _::Correlation{T,U}, chan) where {T,U<:OnePoint{T}} =
    modular_param(chan, x)

conj_q(x, _::Correlation{T,U}) where {T,U<:FourPoints} = conj(x)
conj_q(τ, _::Correlation{T,U}) where {T,U<:OnePoint} = -conj(τ)

total_prefactor(b::ChiralBlock{T,U}, x::PositionCache) where {T,U<:FourPoints} =
    x.prefactor * (x.q_powers[2])^b.channel_dimension.δ
total_prefactor(b::ChiralBlock{T,U}, x::PositionCache) where {T,U<:OnePoint} =
    x.prefactor * x.q^b.channel_dimension.δ
total_prefactor(b, x::Number) = total_prefactor(b, PositionCache(x, b))

function evaluate(b::ChiralBlock{T,U}, x::PositionCache)::T where {T,U}
    d = b.channel_dimension
    p = total_prefactor(b, x)
    h = evaluate_series(b, x)

    if isdegenerate(d)
        # add log(q or 16q) * \sum C^N_rs (q or 16q)^N
        h += x.logq * evalpoly(x, b._missing_terms)
    end

    return p * h
end

function evaluate_der(b::ChiralBlock{T,U}, x::PositionCache)::T where {T,U}
    d = b.channel_dimension
    qor16q = x.q_powers[2]
    p = x.prefactor * (qor16q)^b.channel_dimension.δ
    h = evaluate_series(b, x)
    hprime = evaluate_series_der(b, x)
    h = muladd(h, 2 * d.P * x.logq, hprime) # H_der = 2*P*log(q or 16q)*H + H'
    return p * h
end

@inline evaluate_lr(bs::LeftRight{ChiralBlock}, x::LRPositionCache) =
    evaluate(bs[:left], x.left), evaluate(bs[:right], x.right)
@inline evaluate_lr_der(bs::LeftRight{ChiralBlock}, x::LRPositionCache) =
    evaluate_der(bs[:left], x.left), evaluate_der(bs[:right], x.right)
@inline evaluate_lr(b::BlockFactorized, x) = evaluate_lr(b.chiral_blocks, x)
@inline evaluate_lr(b::LogarithmicBlock, x) = evaluate_lr(b.chiral_blocks, x)
@inline evaluate_lr_op(b::LogarithmicBlock, x) = evaluate_lr(b.chiral_blocks_op, x)
@inline evaluate_lr_der(b::LogarithmicBlock, x) = evaluate_lr_der(b.chiral_blocks_der, x)

@inline function evaluate(b::BlockFactorized{T,U}, x::LRPositionCache)::T where {T,U}
    return prod(evaluate_lr(b, x))
end

function evaluate(b::LogarithmicBlock{T,U}, x::LRPositionCache)::T where {T,U}
    V = b.channel_field
    r, s = V.indices
    s < 0 && return zero(T) # by convention G_(r, s<0) = 0
    Prs = V.P[:left]

    Freg, Fbar = evaluate_lr(b, x)
    F, Fregbar = evaluate_lr_op(b, x)
    R, Rbar = b.R, b.Rbar

    if isaccidentallynonlogarithmic(b)
        Freg * Fbar + (-1)^b.nbzeros * R / Rbar * F * Fregbar

    elseif islogarithmic(b)
        Fder, Fderbar = evaluate_lr_der(b, x)

        l = b.ell

        (Freg - R / 2 / Prs * Fder) * Fbar +
        R / Rbar * F * (Fregbar - Rbar / 2 / Prs * Fderbar) +
        R / 2 / Prs * l * F * Fbar
    end
end

function evaluate(b::InterchiralBlock{T,U}, x)::T where {T,U}
    res = zero(T)
    for i in eachindex(b.blocks)
        res += evaluate((b.blocks)[i], x) .* b.shifts[i]
    end
    return res
end

function evaluate(b::LinearCombinationBlock{T,U}, x)::T where {T,U}
    res = zero(T)
    for i in eachindex(b.blocks)
        res += evaluate((b.blocks)[i], x) .* b.coefficients[i]
    end
    return res
end

evaluate(b::ChiralBlock, x::Number) = evaluate(b, PositionCache(x, b))
evaluate(b::NonChiralBlock, x::Number) = evaluate(b, LRPositionCache(x, b))
evaluate_der(b::ChiralBlock, x::Number) = evaluate_der(b, PositionCache(x, b))
evaluate_der(b::NonChiralBlock, x::Number) = evaluate_der(b, LRPositionCache(x, b))
