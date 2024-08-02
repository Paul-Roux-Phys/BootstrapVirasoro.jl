"""Factor \ell_{(r,s)} that appears in logarithmic blocks"""
function ell(corr::FourPointCorrelation, r, s)
    c = corr.charge
    B, β = c.B, c.β
    βm1P_ext = [[corr.fields[i].P[left]/β for i in 1:4], [corr.fields[i].P[right]/β for i in 1:4]]

    term1(j) = digamma_reg(-2*βm1P(B, r, j)) + digamma_reg(2*βm1P(B, r, -j))

    res = -big(4)*oftype(B, π)/tan(oftype(B, π)*s/B)

    term3(j, lr, pm1, pm2, a, b) = digamma_reg(1/2 + (lr == left ? -1 : 1)*βm1P(B, r, j) +
                                     pm1*βm1P_ext[lr][a] + pm2*βm1P_ext[lr][b])

    return res + 4*sum(term1(j) for j in 1-s:s) -
        sum(term3(j, lr, pm1, pm2, a, b)
                        for pm1 in (-1,1)
                        for pm2 in (-1,1)
                        for j in 1-s:2:s-1
                        for (a,b) in ((1,2), (3, 4))
                        for lr in (left, right)
        )
end

function ell(corr::OnePointCorrelation, r, s)
    println("ell(corr::OnePointCorrelation) not implemented yet")
    return 0
end

function evaluate_logarithmic(
    c::CentralCharge,
    corr::Correlation,
    b::Block,
    pos::Number,
    s::Symbol,
    der,
    reg
)
    b1 = b
    b2 = FourPointBlock(b.Nmax, c, corr, b.channel, swap_lr(b.channel_field))
    # block with momenta (P_(r,-s), P_(r,s)) in the channel

    F_Prms = evaluate_chiral(c, corr, b2, left) # F_{P_(r,-s)}
    F_Prms_bar = evaluate_chiral(conj(x_chan), b1, right) # \bar F_{P_(r,-s)}
    F_der_Prms = evaluate_chiral(x_chan, b2, left, der=true) # F'_{P_(r,-s)}
    F_der_Prms_bar = evaluate_chiral(conj(x_chan), b1, right, der=true) # \bar F'_{P_(r,-s)}
    F_reg_Prs = evaluate_chiral(x_chan, b1, left, reg=true) # F^reg_{P_(r,s)}
    F_reg_Prs_bar = evaluate_chiral(conj(x_chan), b2, right, reg=true) # \bar F^reg_{P_(r,s)}

    R     = corr._Rmn[b.channel][(r, s)][left] # Vchan.P[left] = P_(r,s)
    R_bar = corr._Rmn[b.channel][(r, s)][right]

    term1 = (F_reg_Prs - R * F_der_Prms) * F_Prms_bar
    term2 = R / R_bar * F_Prms * (F_reg_Prs_bar - R_bar * F_der_Prms_bar)
    term3 = -R * ell(corr, r, s) * F_Prms * F_Prms_bar

    # return F_Prms, F_Prms_bar, F_der_Prms, F_der_Prms_bar, F_reg_Prs, F_reg_Prs_bar, ell(corr, r, s), R, R_bar
    return channel_sign(block, x) * (term1 + term2 + term3)
end

