series_argument(x, V::Union{FourFields, FourDimensions}) = 16*qfromx(x)
series_argument(τ, V::Union{OneField, OneDimension}) = exp(im*oftype(τ, π)*τ)
series_argument(x, b::BlockChiral) = series_argument(x, b.corr.dims)
series_argument(x, b::BlockNonChiral) = series_argument(x, b.corr.fields)

function evaluate_series(b::Block, q; der=false)
    der && return evalpoly(q, b._coefficients_der)
    evalpoly(q, b._coefficients)
end

get_position(x, V::Union{FourFields, FourDimensions}, b::Block) = crossratio(b.channel, x)
get_position(x, V::Union{OneField, OneDimension}, b::Block) = x
get_position(x, c::CorrelationChiral, b::Block) = get_position(x, c.dims, b)
get_position(x, c::Correlation, b::Block) = get_position(x, c.fields, b)

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
        missingterm = [
            (N, r, s) in keys(b._CNmn) ? b._CNmn[(N, r, s)] : zero(x)
            for N in 0:c.Nmax
        ]
        h += log(q) * evalpoly(q, missingterm)
    end

    p * h
end

function evaluate(b::BlockNonChiral, x, lr; der=false)
    xs = (x, conj(x))
    if der
        return evaluate(b.chiral_blocks_der[lr], xs[lr], der=der)
    end
    evaluate(b.chiral_blocks[lr], xs[lr])
end

function evaluate(b::BlockNonChiral, x)
    if !(isaccidentallynonlogarithmic(b) || islogarithmic(b))
        return prod(evaluate(b, x, lr) for lr in (:left, :right))
    else
        V      = b.channel_field
        r, s   = V.indices
        Pp = ConformalDimension(V.c, r=r, s=s).P
        βm1Prs = (r+s/V.c.B)/2
        b_op   = Block(b.corr, b.channel, swap_lr(b.channel_field), b.Nmax)
        if s < 0
            b_op, b = b, b_op # b has left, right dims (P_(r, s>0), P_(r, -s<0))
                              # b_op  has dims (P_(r, -s<0), P_(r, s>0))
        end

        Freg    = evaluate(b, x, :left)
        Fbar    = evaluate(b, x, :right)
        F       = evaluate(b_op, x, :left)
        Fregbar = evaluate(b_op, x, :right)

        if isaccidentallynonlogarithmic(b)
            Rreg = b.corr._Rmn_reg[:left][b.channel][(r, s)]
            Rregbar = b.corr._Rmn_reg[:right][b.channel][(r, s)]
            nbzeros = Rmn_zero_order(r, s, b.chiral_blocks[:left].dims)
            return Freg*Fbar + (-1)^nbzeros * Rreg/Rregbar * F*Fregbar

        elseif islogarithmic(b)
            Fder = evaluate(b_op, x, :left, der=true)
            Fderbar = evaluate(b, x, :right, der=true)

            R = b.corr._Rmn[:left][b.channel][(r, s)]
            Rbar = b.corr._Rmn[:right][b.channel][(r, s)]

            return (
                (Freg - R / 2 / Pp * Fder) * Fbar +
                R / Rbar * F * (Fregbar - Rbar / 2 / Pp * Fderbar) -
                R * ell(b, r, s) / V.c.B / 2 / βm1Prs * F * Fbar
            )
        end
    end
end
