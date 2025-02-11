"""Factor ``\\ell_{(r,s)}`` that appears in logarithmic blocks"""
ell(b::Block, r, s) = ell(b.corr.fields, b.channel, r, s)

function ell(V::FourFields, chan, r, s)
    V = permute_fields(V, chan)
    c = V[1].c
    βm1P_ext = Tuple(
        Tuple(V[i].P[lr] / c.β for i in 1:4)
        for lr in (:left, :right)
    )
    res = -4 * oftype(c.B, π) / tan(oftype(c.B, π) * s / c.B)
    res += 4 * sum(
        digamma_reg(-2 * Prs(r, j, c) / c.β) + digamma_reg(2 * Prs(r, -j, c) / c.β)
        for j in 1-s:s
    )
    res -= sum(
        digamma_reg(1 // 2 + (lr == :left ? -1 : 1) * Prs(r, j, c)/c.β +
                    pm1 * βm1P_ext[lr][a] + pm2 * βm1P_ext[lr][b])
        for pm1 in (-1, 1)
        for pm2 in (-1, 1)
        for j in 1-s:2:s-1
        for (a, b) in ((1, 2), (3, 4))
        for lr in (:left, :right)
    )
    return res
end

function ell(V::OneField, chan, r, s)
    V1 = V[1]
    c = V1.c
    res = -4 * oftype(c.B, π) / tan(oftype(c.B, π) * s / c.B)
    res += 4 * sum(
        digamma_reg(-2 * Prs(r, j, c) / c.β) + digamma_reg(2 * Prs(r, -j, c) / c.β)
        for j in 1-s:s
    )
    res += 2 * sum(
        digamma_reg((Prs(1, pm1, c) + pm2 * V[1].P[:right] + 2 * Prs(r, j, c)) / c.β) +
        digamma_reg((Prs(1, pm1, c) + pm2 * V[1].P[:left] - 2 * Prs(r, j, c)) / c.β)
        for j in 1-s:2:s-1
        for pm1 in (-1, 1)
        for pm2 in (-1, 1)
    )
    return res
end

function isaccidentallynonlogarithmic(b::BlockNonChiral)
    V = b.channel_field
    !(V.isKac && V.r%1 == V.s%1 == 0) && return false
    r, s = V.indices
    return Rmn_zero_order(r, s, permute_dimensions(b.chiral_blocks[:left].dims, b.channel)) > 0
end

function islogarithmic(b::BlockNonChiral)
    V = b.channel_field;
    V.isKac && V.r%1 == V.s%1 == 0 && V.r*V.s != 0 && !(isaccidentallynonlogarithmic(b))
end
