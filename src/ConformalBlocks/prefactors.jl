#===========================================================================================
Sphere prefactor
===========================================================================================#

"""Nome `q` from the cross-ratio `x`"""
@memoize qfromx(x) = exp(- (π * ellipticK(1 - x) / ellipticK(x)))

"""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""
    channelprefactor_chiral(d::FourDimensions, b, x)

Prefactor to get t- or u-channel blocks from the s-channel block
"""
function channelprefactor_chiral(d::FourDimensions, chan, x)
    chan === :u && return x^(2*d[1].Δ)
    return 1
end

"""Sign (-1)^{S_1+S_2+S_3+S_4} when changing from s to t or u channels"""
function channel_sign(b::Block, x)
    b.channel === :s && return 1
    return 1 # (-1)^(sum(spin.(corr.fields)))
end

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
@memoize function crossratio(chan, x)
    chan === :s && return x
    chan === :t && return 1 .- x
    chan === :u && return 1 ./ x
    error(
        """Incorrect channel specification in crossratio(channel, x):
        must be in $channels"""
    )
end

"""
    blockprefactor_chiral(d::FourDimensions, b::BlockChiral, x)

Prefactor for getting the chiral block F from H. 
The argument `lr` indicates if we are working with a left or right moving block
"""
@memoize function prefactor_chiral(d::FourDimensions, chan::Symbol, x)
    ds = permute_dimensions(d, chan)
    c = d[1].c.c

    a = (c-1)/24
    e0 = -ds[1].δ - ds[2].δ - a
    e1 = -ds[1].δ - ds[4].δ - a
    e2 = sum(ds[i].δ for i in 1:4) + a

    q = qfromx(x)
    
    return channelprefactor_chiral(d, chan, x) *
        x^e0 * (1 - x)^e1 * jtheta3(0, q)^(-4 * e2), q
end

function blockprefactor_chiral(d::FourDimensions, b, x)
    q = qfromx(x)
    pref = prefactor_chiral(d, b.channel, x)
    q = pref[2]
    pref[1] * (16 * q)^b.channel_dimension.δ
end


#===========================================================================================
Torus prefactor
===========================================================================================#
@memoize qfromτ(τ) = exp(im*(π*τ))

function blockprefactor_chiral(d::OneDimension, b::BlockChiral, τ)
    q = qfromτ(τ)
    q^b.channel_dimension.δ / etaDedekind(complex(τ))
end

blockprefactor_chiral(b::BlockChiral, x) = blockprefactor_chiral(b.corr.dims, b, x)
