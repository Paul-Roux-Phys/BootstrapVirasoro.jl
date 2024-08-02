"""Nome `q` from the cross-ratio `x`"""
qfromx(x) = exp(-oftype(x, π) * ellipticK(1 - x) / ellipticK(x))

""""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#
"""Prefactor to get t- or u-channel blocks from the s-channel block"""
function channelprefactor_chiral(
    corr::FourPointCorrelation{T},
    block::FourPointBlock{T},
    x,
    lr
) where {T}
    block.channel === :u && return x^(2*corr.fields[1].Δ[lr])
    return one(T)
end

"""Sign (-1)^{S_1+S_2+S_3+S_4} when changing from s to t or u channels"""
function channel_sign(block::FourPointBlock, x)
    channel = block.channel
    channel === :s && return 1
    return 1 # (-1)^(sum(spin.(corr.fields)))
end

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
function crossratio(channel, x)
    channel === :s && return x
    channel === :t && return 1-x
    channel === :u && return 1/x
    error(
        """Incorrect channel specification in crossratio(channel, x):
        must be in $channels""")
end

#===========================================================================================
Block prefactor
===========================================================================================#
"""
Prefactor for getting the chiral block F from H. 
The argument `lr` indicates if we are working with a left or right moving block
"""
function blockprefactor_chiral(
    cc::CentralCharge,
    corr::FourPointCorrelation,
    b::FourPointBlock,
    x,
    lr)

    Vs = permute_fields(corr.fields, b.channel)

    c = cc.c
    e0 = -Vs[1].δ[lr] - Vs[2].δ[lr] - (c - 1) / 24
    e1 = -Vs[1].δ[lr] - Vs[4].δ[lr] - (c - 1) / 24
    e2 = sum(Vs[i].δ[lr] for i in 1:4) + (c - 1) / 24
    q = qfromx(x)

    res = x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2) * (16 * q)^b.channel_field.δ[lr]

    return channelprefactor_chiral(corr, b, x, lr) * res
end