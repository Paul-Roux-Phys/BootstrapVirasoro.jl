qfromx(x) = exp(- (π * ellipticK(1 - x) / ellipticK(x)))
xfromq(q) = jtheta2(0, q)^4 / jtheta3(0, q)^4
qfromτ(τ) = exp(2*im*(π*τ))
τfromx(x) = (log(qfromx(x)) / π) / im

struct PosCache{T}
    x::T
    prefactor::T
    q::T
    logq::T
    q_powers::Vector{T}
end

const LRPosCache = LeftRight{PosCache}

function PosCache(x, ds::NTuple{4,CD{T}}, chan::Symbol, Δmax) where {T}
    ds = permute_4(ds, chan)

    q = qfromx(x)

    e0 = -ds[1].Δ - ds[2].δ
    chan === :u && (e0 += 2ds[1].Δ)
    e1 = -ds[1].Δ - ds[4].δ
    e2 = sum(ds[i].δ for i = 1:3) + ds[4].Δ

    prefactor = x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2)

    sq = 16q
    q_powers = ones(T, Δmax+1)
    for i = 2:(Δmax+1)
        q_powers[i] = q_powers[i-1] * sq
    end

    return PosCache{T}(x, prefactor, q, log(sq), q_powers)
end

function PosCache(τ, _::NTuple{1,CD{T}}, _::Symbol, Δmax) where {T}
    q = qfromτ(τ)
    prefactor = 1 / etaDedekind(complex(τ))
    q_powers = ones(T, Δmax+1)
    for i = 2:(Δmax+1)
        q_powers[i] = q_powers[i-1] * q
    end

    return PosCache{T}(τ, prefactor, q, log(q), q_powers)
end

PosCache(x, co::CCo, chan) = PosCache(x, co.fields, chan, co.Δmax)
PosCache(x, b::CBlock) = PosCache(x, b.corr, b.chan)

function LeftRight{PosCache}(x, co::NCCo{T}, chan) where {T}
    xbar = conj_q(x, co)
    return LRPosCache(PosCache(x, co[:left], chan), PosCache(xbar, co[:right], chan))
end
LeftRight{PosCache}(x, b::NCBlock) = LRPosCache(x, b.corr, b.chan)

function evalpoly(x::PosCache, coeffs::Vector{T}) where {T}
    res = zero(T)
    for i = 1:length(coeffs)
        res += coeffs[i] * x.q_powers[i]
    end
    return res
end

eval_series(b::CBlock, x::PosCache) = evalpoly(x, b.coeffs)
eval_series_der(b::CBlock, x::PosCache) = evalpoly(x, b.coeffs_der)
eval_series(b::CBlock, x::Number) = eval_series(b, PosCache(x, b.corr, b.chan))
eval_series_der(b::CBlock, x::Number) =
    eval_series_der(b, PosCache(x, b.corr, b.chan))

conj_q(x, _::Correlation4) = conj(x)
conj_q(τ, _::Correlation1) = -conj(τ)

prefactor(b::ChiralBlock, x::Number) = PosCache(x, b).prefactor
prefactor(b::NonChiralBlock, x::LRPosCache) =
    prefactor(b.cblocks.left, x.left) * prefactor(b.cblocks.right, x.right)
prefactor(b::NonChiralBlock, x::Number) = prefactor(b, LRPosCache(x, b))

total_prefactor(b::CBlock, x::PosCache, _::Correlation4) =
    x.prefactor * (x.q_powers[2])^b.chan_dim.δ
total_prefactor(b::CBlock, x::PosCache, _::Correlation1) =
    x.prefactor * x.q^b.chan_dim.δ
total_prefactor(b::CBlock, x::PosCache) = total_prefactor(b, x, b.corr)
total_prefactor(b::CBlock, x::Number) = total_prefactor(b, PosCache(x, b))

function (b::CBlock{T})(x::PosCache)::T where {T}
    d = b.chan_dim
    p = total_prefactor(b, x)
    h = eval_series(b, x)

    if d.degenerate
        # add log(q or 16q) * \sum C^N_rs (q or 16q)^N
        h += x.logq * evalpoly(x, b.missing_terms)
    end

    return p * h
end

function eval_der(b::CBlock{T}, x::PosCache)::T where {T}
    d = b.chan_dim
    qor16q = x.q_powers[2]
    p = x.prefactor * (qor16q)^b.chan_dim.δ
    h = eval_series(b, x)
    hprime = eval_series_der(b, x)
    h = muladd(h, 2 * d.P * x.logq, hprime) # H_der = 2*P*log(q or 16q)*H + H'
    return p * h
end

eval_lr(bs::LR{CBlock{T}}, x) where {T} =
    bs.left(x.left), bs.right(x.right)
eval_lr_der(bs::LeftRight{CBlock{T}}, x) where {T} =
    eval_der(bs.left, x.left), eval_der(bs.right, x.right)
eval_lr(b::FactorizedBlock{T}, x) where {T} = eval_lr(b.cblocks, x)
eval_lr(b::LogBlock, x) = eval_lr(b.cblocks, x)
eval_lr_op(b::LogBlock, x) = eval_lr(b.cblocks_op, x)
eval_lr_der(b::LogBlock, x) = eval_lr_der(b.cblocks_der, x)

function (b::FactorizedBlock{T})(x::LRPosCache)::T where {T}
    lr = eval_lr(b, x)
    return lr[1] * lr[2]
end

function (b::LogBlock{T})(x::LRPosCache)::T where {T}
    V = b.chan_field
    _, s = indices(V)
    s < 0 && return zero(T) # by convention G_(r, s<0) = 0
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

function (b::LCBlock{T})(x)::T where {T}
    res = zero(T)
    for i in eachindex(b.blocks)
        res += (b.blocks)[i](x) .* b.coeffs[i]
    end
    return res
end

(b::CBlock)(x::Number) = b(PosCache(x, b))
(b::NCBlock)(x::Number) = b(LRPosCache(x, b))
eval_der(b::CBlock, x::Number) = eval_der(b, PosCache(x, b))
eval_der(b::NCBlock, x::Number) = eval_der(b, LRPosCache(x, b))
