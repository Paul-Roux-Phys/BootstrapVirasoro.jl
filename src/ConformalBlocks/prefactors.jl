#===========================================================================================
Sphere prefactor
===========================================================================================#

"""Nome `q` from the cross-ratio `x`"""
@inline qfromx(x) = exp(- (π * ellipticK(1 - x) / ellipticK(x)))

"""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""

function channelprefactor_chiral(d::FourDimensions, chan, x)
    chan === :u && return x^(2*d[1].Δ)
    return 1
end

"""
    blockprefactor_chiral(d::FourDimensions, b::BlockChiral, x)

Prefactor for getting the chiral block F from H. 
The argument `lr` indicates if we are working with a left or right moving block
"""
function prefactor_chiral(d::FourDimensions, chan::Symbol, x)
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
    pref, q = prefactor_chiral(d, b.channel, x)
    pref * (16 * q)^b.channel_dimension.δ
end


#===========================================================================================
Torus prefactor
===========================================================================================#
qfromτ(τ) = exp(im*(π*τ))

function prefactor_chiral(d::OneDimension, chan::Symbol, τ)
    q = qfromτ(τ)
    1 / etaDedekind(complex(τ)), q
end

function blockprefactor_chiral(d::OneDimension, b, τ)
    pref, q = prefactor_chiral(d, b.channel, τ)
    pref * q^b.channel_dimension.δ 
end

blockprefactor_chiral(b::BlockChiral, x) = blockprefactor_chiral(b.corr.dims, b, x)
