βm1P(B, r, s) = 1/2*(r+s/B) # \beta^{-1}P_{(r,s)}

"""Factor \ell_{(r,s)} that appears in logarithmic blocks"""
function ell(V::FourFields, chan, r, s)
    V = permute_fields(V, chan)
    B, β = V[1].c.B, V[1].c.β
    βm1P_ext = Tuple(
        Tuple(V[i].P[lr]/β for i in 1:4)
        for lr in (:left, :right)
    )
    res = -big(4)*oftype(B, π)/tan(oftype(B, π)*s/B)
    return res +
        4*sum(
            digamma_reg(-2*βm1P(B, r, j)) + digamma_reg(2*βm1P(B, r, -j))
            for j in 1-s:s
        ) - sum(
            digamma_reg(1//2 + (lr == :left ? -1 : 1)*βm1P(B, r, j) +
                                                  pm1*βm1P_ext[lr][a] + pm2*βm1P_ext[lr][b])
            for pm1 in (-1,1)
            for pm2 in (-1,1)
            for j in 1-s:2:s-1
            for (a,b) in ((1,2), (3, 4))
            for lr in (:left, :right)
        )
end

function ell(V::OneField, channel, r, s)
    println("ell(V::OneField, r, s) not implemented yet")
    return 0
end

ell(b::Block, r, s) = ell(b.corr.fields, b.channel, r, s)

function isaccidentallynonlogarithmic(b::BlockNonChiral)
    V = b.channel_field
    !(V.isKac && V.r%1 == V.s%1 == 0) && return false
    r, s = V.indices
    return false
    # return (r, s) in keys(b.corr._Rmn[b.channel][:left]) ||
    #        (r, s) in keys(b.corr._Rmn[b.channel][:right])
end

function islogarithmic(b::BlockNonChiral)
    V = b.channel_field;
    V.isKac && V.r%1 == V.s%1 == 0 && V.r*V.s != 0
end