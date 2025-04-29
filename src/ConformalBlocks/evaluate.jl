series_argument(x, V::Union{FourFields, FourDimensions}) = 16*qfromx(x)
series_argument(τ, V::Union{OneField, OneDimension}) = exp(im*oftype(τ, π)*τ)
series_argument(x, b::BlockChiral) = series_argument(x, b.corr.dims)
series_argument(x, b::BlockNonChiral) = series_argument(x, b.corr.fields)

function evaluate_series(b::BlockChiral{T}, q; der=false) where {T}
    der && return evalpoly_buf(q, b._coefficients_der)
    evalpoly_buf(q, b._coefficients)
end

get_position(x, V::Union{FourFields, FourDimensions}, b::Block) = crossratio(b.channel, x)
get_position(x, V::Union{OneField, OneDimension}, b::Block) = x
get_position(x, c::CorrelationChiral, b::Block) = get_position(x, c.dims, b)
get_position(x, c::Correlation, b::Block) = get_position(x, c.fields, b)

conj_q(x, V::Union{FourFields, FourDimensions}, b::Block) = conj.(x)
conj_q(τ, V::Union{OneField, OneDimension}, b::Block) = @. -conj(τ)
conj_q(x, c::CorrelationChiral, b::Block) = conj_q.(x, Ref(c.dims), Ref(b))
@memoize conj_q(x, b::Block) = conj_q(x, b.corr[:left], b)

function evaluate(b::BlockChiral, x::T; der=false)::T where {T}
    co = b.corr
    y = get_position(x, co.dims, b)
    q = series_argument.(y, Ref(b))
    
    d = b.channel_dimension
    p = blockprefactor_chiral.(Ref(co.dims), Ref(b), y)

    h = evaluate_series(b, q)

    # add the q-dependent parts
    if der
        hprime = evaluate_series(b, q, der=true)
        h = muladd.(h, 2*d.P*log.(q), hprime) # H_der = 2*P*log(q or 16q)*H + H'
    elseif d.isKac && d.r%1 == d.s%1 == 0 && d.r > 0 && d.s > 0
        r, s = d.indices
        # add log(q or 16q) * \sum C^N_rs (q or 16q)^N
        missingterm = [
            (N, r, s) in keys(b._CNmn) ? b._CNmn[(N, r, s)] : zero(eltype(values(b._CNmn)))
            for N in 0:co.Nmax
        ]
        h += log.(q) .* evalpoly_buf(q, missingterm)
    end

    p .* h
end

function evaluate(b::BlockFactorized, x::T, lr)::T where {T}
    xs = (x, conj_q(x, b))
    evaluate(b.chiral_blocks[lr], xs[lr])
end

function evaluate(b::BlockLogarithmic, x::T, lr; der=false, op=false)::T where {T}
    xs = (x, conj_q(x, b))
    der && return evaluate(b.chiral_blocks_der[lr], xs[lr], der=der)
    op && return evaluate(b.chiral_blocks_op[lr], xs[lr])
    evaluate(b.chiral_blocks[lr], xs[lr])
end

function evaluate(b::BlockFactorized, x::T)::T where {T}
    evaluate(b[:left], x) .* evaluate(b[:right], x)
end

function evaluate(b::BlockLogarithmic, x::T)::T where {T}
    V = b.channel_field
    r, s = V.indices
    s < 0 && return x isa AbstractArray ? zeros(eltype(x), size(x)) : zero(x) # by convention G_(r, s<0) = 0
    Prs = V.P[:left]

    Freg = evaluate(b[:left], x)
    Fbar = evaluate(b[:right], x)
    F = evaluate(b, x, :left, op=true)
    Fregbar = evaluate(b, x, :right, op=true)

    if isaccidentallynonlogarithmic(b)
        Rreg = b.corr._Rmn_reg[:left][b.channel][(r, s)]
        Rregbar = b.corr._Rmn_reg[:right][b.channel][(r, s)]
        nbzeros = Rmn_zero_order(r, s, b.chiral_blocks[:left].dims)
        @. Freg * Fbar + (-1)^nbzeros * Rreg / Rregbar * F * Fregbar

    elseif islogarithmic(b)
        Fder = evaluate(b, x, :left, der=true)
        Fderbar = evaluate(b, x, :right, der=true)

        R = b.corr._Rmn[:left][b.channel][(r, s)]
        Rbar = b.corr._Rmn[:right][b.channel][(r, s)]
        l = ell(b)

        @. (Freg - R / 2 / Prs * Fder) * Fbar +
        R / Rbar * F * (Fregbar - Rbar / 2 / Prs * Fderbar) +
        R / 2 / Prs * l * F * Fbar
    end
end

function evaluate(b::BlockInterchiral, x::T)::T where {T}
    sum(
        b.shifts[i] .* evaluate((b.blocks)[i], x)
        for i in eachindex(b.blocks)
    )
end
