"""Factor ``\\ell_{(r,s)}`` that appears in logarithmic blocks"""
ell(b::Block, r, s) = ell(b.corr.fields, b.channel, r, s)

"""Factor \ell_{(r,s)} that appears in logarithmic blocks"""
function ell(V::FourFields, chan, r, s)
    V = permute_fields(V, chan)
    β = V[1].c.β
    res = 4 * (π / tan( s * (π / β^2)))
    res += 4 * sum(
        digamma_reg(-2 * Prs(r, j, β) / β) + digamma_reg(2 * Prs(r, -j, β) / β)
        for j in 1-s:s
    )
    res -= sum(
        digamma_reg(1 // 2 + (lr == :left ? -1 : 1) *
            (Prs(r, j, β) + pm1 * V[a].P[lr] + pm2 * V[b].P[lr]) / β
        )
        for pm1 in (-1, 1)
        for pm2 in (-1, 1)
        for j in 1-s:2:s-1
        for (a, b) in ((1, 2), (3, 4))
        for lr in (:left, :right)
    )
    res / β
end
    
function ell(V::OneField, chan, r, s)
    β = V[1].c.β
    res = 4 * (π / tan( s * (π / β^2)))
    res += 4 * sum(
        digamma_reg(-2 * Prs(r, j, β) / β) + digamma_reg(2 * Prs(r, -j, β) / β)
        for j in 1-s:s
    )
    res -= 2 * sum(
        digamma_reg(1//2 + pm1 * inv(β)^2 / 2 +
            (pm2 * V[1].P[lr] + 2 * (lr == :left ? -1 : 1) * Prs(r, j, β)) / β)
        for j in 1-s:2:s-1
        for pm1 in (-1, 1)
        for pm2 in (-1, 1)
        for lr in (:left, :right)
    )
    res / β
end

ell(c::Correlation, chan, r, s) = ell(c.fields, chan, r, s)
ell(b::Block) = ell(b.corr, b.channel, b.channel_field.r, b.channel_field.s)

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
