"""Nome `q` from the cross-ratio `x`"""
qfromx(x) = exp(-oftype(x, π) * ellipticK(1 - x) / ellipticK(x))

""""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#
"""Prefactor to get t- or u-channel blocks from the s-channel block"""
function channelprefactor_chiral(V::FourFields, b, x, lr)
    b.channel === :u && return x^(2*V[1].Δ[lr])
    return 1
end

"""Sign (-1)^{S_1+S_2+S_3+S_4} when changing from s to t or u channels"""
function channel_sign(b::Block, x)
    b.channel === :s && return 1
    return 1 # (-1)^(sum(spin.(corr.fields)))
end

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
function crossratio(chan, x)
    chan === :s && return x
    chan === :t && return 1-x
    chan === :u && return 1/x
    error(
        """Incorrect channel specification in crossratio(channel, x):
        must be in $channels"""
    )
end

#===========================================================================================
Block prefactor
===========================================================================================#
"""
Prefactor for getting the chiral block F from H. 
The argument `lr` indicates if we are working with a left or right moving block
"""
function blockprefactor_chiral(V::FourFields, b::BlockChiral, x)
    Vs = permute_fields(V, b.channel)
    lr = b.lr

    a = (V[1].c.c-1)/24
    e0 = -Vs[1].δ[lr] - Vs[2].δ[lr] - a
    e1 = -Vs[1].δ[lr] - Vs[4].δ[lr] - a
    e2 = sum(Vs[i].δ[lr] for i in 1:4) + a
    q = qfromx(x)
    
    return channelprefactor_chiral(V, b, x, lr) *
           x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2) * (16 * q)^b.channel_dimension.δ
end
