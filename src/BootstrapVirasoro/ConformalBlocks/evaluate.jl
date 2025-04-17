series_argument(x, V::Union{FourFields, FourDimensions}) = 16*qfromx(x)
series_argument(τ, V::Union{OneField, OneDimension}) = exp(im*oftype(τ, π)*τ)
series_argument(x, b::BlockChiral) = series_argument(x, b.corr.dims)
series_argument(x, b::BlockNonChiral) = series_argument(x, b.corr.fields)

function evaluate_series(b::BlockChiral, q; der=false)
    der && return evalpoly(q, b._coefficients_der)
    evalpoly(q, b._coefficients)
end

get_position(x, V::Union{FourFields, FourDimensions}, b::Block) = crossratio(b.channel, x)
get_position(x, V::Union{OneField, OneDimension}, b::Block) = x
get_position(x, c::CorrelationChiral, b::Block) = get_position(x, c.dims, b)
get_position(x, c::Correlation, b::Block) = get_position(x, c.fields, b)

conj_q(x, V::Union{FourFields, FourDimensions}, b::Block) = conj(x)
conj_q(τ, V::Union{OneField, OneDimension}, b::Block) = complex(-real(τ), imag(τ))
conj_q(x, c::CorrelationChiral, b::Block) = conj_q(x, c.dims, b)
conj_q(x, b::Block) = conj_q(x, b.corr[:left], b)

function evaluate(b::BlockChiral, x; der=false)
    c = b.corr
    y = get_position(x, c.dims, b)
    q = series_argument(y, b)
    
    d = b.channel_dimension
    p = blockprefactor_chiral(c.dims, b, y)

    h = evaluate_series(b, q)

    # add the q-dependent parts
    if der
        hprime = evaluate_series(b, q, der=true)
        h = muladd(h, 2*d.P*log(q), hprime) # H_der = 2*P*log(q or 16q)*H + H'
    elseif d.isKac && d.r%1 == d.s%1 == 0 && d.r > 0 && d.s > 0
        r, s = d.indices
        # add log(q or 16q) * \sum C^N_rs (q or 16q)^N
        missingterm = [
            (N, r, s) in keys(b._CNmn) ? b._CNmn[(N, r, s)] : zero(x)
            for N in 0:c.Nmax
        ]
        h += log(q) * evalpoly(q, missingterm)
    end

    p * h
end

function evaluate(b::BlockFactorized, x, lr)
    xs = (x, conj_q(x, b))
    evaluate(b.chiral_blocks[lr], xs[lr])
end

function evaluate(b::BlockLogarithmic, x, lr; der=false, op=false)
    xs = (x, conj_q(x, b))
    der && return evaluate(b.chiral_blocks_der[lr], xs[lr], der=der)
    op && return evaluate(b.chiral_blocks_op[lr], xs[lr])
    evaluate(b.chiral_blocks[lr], xs[lr])
end

function evaluate(b::BlockFactorized, x)
    prod(evaluate(b, x, lr) for lr in (:left, :right))
end

function evaluate(b::BlockLogarithmic, x; debug=false)
    V = b.channel_field
    r, s = V.indices
    s < 0 && return 0 # by convention G_(r, s<0) = 0
    P_rs = Prs(r, s, V.c)

    Freg = evaluate(b, x, :left)
    Fbar = evaluate(b, x, :right)
    F = evaluate(b, x, :left, op=true)
    Fregbar = evaluate(b, x, :right, op=true)

    if isaccidentallynonlogarithmic(b)
        Rreg = b.corr._Rmn_reg[:left][b.channel][(r, s)]
        Rregbar = b.corr._Rmn_reg[:right][b.channel][(r, s)]
        nbzeros = Rmn_zero_order(r, s, b.chiral_blocks[:left].dims)
        return Freg * Fbar + (-1)^nbzeros * Rreg / Rregbar * F * Fregbar

    elseif islogarithmic(b)
        Fder = evaluate(b, x, :left, der=true)
        Fderbar = evaluate(b, x, :right, der=true)

        R = b.corr._Rmn[:left][b.channel][(r, s)]
        Rbar = b.corr._Rmn[:right][b.channel][(r, s)]

        terms = (
            (Freg - R / 2 / P_rs * Fder) * Fbar,
            R / Rbar * F * (Fregbar - Rbar / 2 / P_rs * Fderbar),
            +R / 2 / P_rs * ell(b) * F * Fbar
        )

        return if debug terms else sum(terms) end
    end
end

function evaluate(b::BlockInterchiral, x)
    sum(
        b.shifts[i] * evaluate((b.blocks)[i], x)
        for i in 1:length(b.blocks)
    )
end
