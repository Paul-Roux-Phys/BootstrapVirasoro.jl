#===========================================================================================

ConformalBlocks.jl contains modules that compute Virasoro four-point conformal blocks on the
sphere and Virasoro one-point conformal blocks on the torus.

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#


"""
Computation of four-point blocks on the sphere.
"""
module FourPointBlocksSphere

export FourPointBlockSphere, block_chiral, block_non_chiral

using ..CFTData, ..FourPointCorrelationFunctions
using Match, EllipticFunctions, Memoization
import ..FourPointCorrelationFunctions: permute_ext_fields, Rmn
import ..JuliVirBootstrap.SpecialFunctions: digamma_reg

#===========================================================================================
Struct FourPointBlockSphere
===========================================================================================#
"""
    FourPointBlockSphere{T}

Composite type that represents the list of arguments of a four-point conformal block:
a channel and a field propagating in the channel. The external fields and central charge are
provided in a `FourPointCorrelation` object.

# Example

```julia-repl
julia> c = CentralCharge("c",0.5); V = Field(c, "δ", 0.6, diagonal = true);
julia> FourPointBlockSphere("s", V)
Four-point block
Channel:        s
Channel Field:
Diagonal field of dimension:
  Δ = 0.5791666666666667 + 0.0im
  P = 0.0 + 0.7745966692414834im
  δ = 0.6000000000000001 + 0.0im
  p = 0.7745966692414834 + 0.0im
```
"""
struct FourPointBlockSphere{T}

    channel::String
    channelField::Field{T}

end

"""Display blocks"""
function Base.show(io::IO, block::FourPointBlockSphere)
    println("Four-point block")
    println("Channel:\t$(block.channel)")
    println("Channel Field:")
    show(block.channelField)
    # println("External Fields:")
    # print("1. "); show(block.extFields[1])
    # print("2. "); show(block.extFields[2])
    # print("3. "); show(block.extFields[3])
    # print("4. "); show(block.extFields[4])
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2

#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#
"""Prefactor to get t- or u-channel blocks from the s-channel block"""
function channelprefactor_chiral(block::FourPointBlockSphere, corr::FourPointCorrelation, x)
    @match block.channel begin
        "s" => 1
        "t" => 1
        "u" => 1/x^(2*corr.fields[1]["Δ"][left])
    end
end

"""Sign (-1)^{S_1+S_2+S_3+S_4} when changing from s to t or u channels"""
function channel_sign(block::FourPointBlockSphere, corr::FourPointCorrelation, x)
    @match block.channel begin
        "s" => 1
        "t" => 1 # (-1)^(sum(spin.(corr.fields)))
        "u" => 1 # (-1)^(sum(spin.(corr.fields)))
    end
end

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
function crossratio(channel, x)
    @match channel begin
        "s" => x
        "t" => 1-x
        "u" => 1/x
    end
end

#===========================================================================================
Set prefactors, relate the cross-ratio x and the elliptic nome q
===========================================================================================#
"""Nome `q` from the cross-ratio `x`"""
qfromx(x) = exp(-π*ellipticK(1-x) / ellipticK(x))

"""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""Prefactor for getting the block F from H. The argument `lr` indicates if we are working
with a left or right moving block"""
function blockprefactor(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

    c = corr.charge["c"]
    e0 = - corr.fields[1]["δ"][lr] - corr.fields[2]["δ"][lr] - (c-1)/24
    e1 = - corr.fields[1]["δ"][lr] - corr.fields[4]["δ"][lr] - (c-1)/24
    e2 = sum(corr.fields[i]["δ"][lr] for i in 1:4) + (c-1)/24
    q=qfromx(x)

    return x^e0 * (1-x)^e1 * jtheta3(0,q)^(-4*e2) * (16*q)^block.channelField["δ"][1]
end

"""Degenerate dimensions"""
δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)

βm1P(B, r, s) = 1/2*(r+s/B) # \beta^{-1}P_{(r,s)}

"""Factor \ell_{(r,s)} that appears in logarithmic blocks"""
function ell(corr, r, s)
    c = corr.charge
    B, β = c["B"], c["β"]
    βm1P_ext = [[corr.fields[i]["P"][left]/β for i in 1:4], [corr.fields[i]["P"][right]/β for i in 1:4]]

    term1(j) = digamma_reg(-2*βm1P(B, r, j)) + digamma_reg(2*βm1P(B, r, -j))

    res = -big(4)*π/tan(π*big(s)/B) # I put big(n)*\pi otherwise n*\pi where n is an integer has double precision instead of bigfloat

    term3(j, lr, pm1, pm2, a, b) = digamma_reg(1/2 + (lr == left ? -1 : 1)*βm1P(B, r, j) + pm1*βm1P_ext[lr][a] + pm2*βm1P_ext[lr][b])

    return res + 4*sum(term1(j) for j in 1-s:s) -
        sum(term3(j, lr, pm1, pm2, a, b)
                        for pm1 in (-1,1)
                        for pm2 in (-1,1)
                        for j in 1-s:2:s-1
                        for (a,b) in ((1,2), (3, 4))
                        for lr in (left, right)
        )
end

#===========================================================================================
Compute the conformal block
===========================================================================================#
function P_squared_ratio_reg(q, c::CentralCharge, V::Field, m, n, lr)
    # check V has integer Kac indices
    if V.isKac && V.r%1 == 0 && V.s%1 == 0 && V.r > 0 && (lr == left && V.s > 0 || lr == right && V.s < 0)
        # if s < 0 and we're computing a right-handed block (\bar F) then the right dimension is P_(r,-s>0)
        P = V["P"][left]
        β = c["β"]
        Pmn = 1/2*(β*m - 1/β*n)
        if V.r == m && V.s == n
            return log(16*q) - 1/(4*P^2)
        else
            return 1/(P^2-Pmn^2)
        end
    else
        error("Trying to compute a regularised block for a field with r=$(V.r) and s=$(V.s) . Both should be positive integers")
    end
end

function block_recursion_coeff(q, c, V, m, n, der, reg, lr)
    β = c["β"]
    P = V["P"][lr]
    Pmn = 1/2*(β*m - 1/β*n)
    if der
        return 2*P*(log(16*q)/(P^2-Pmn^2) - 2*P/(P^2-Pmn^2)^2) # 2P (log(16q)/(δ-\delta_{m,n}) - 1/(δ-\delta_{m,n})^2)
    elseif reg
        return P_squared_ratio_reg(q, c, V, m, n, lr) # log(16q) - 1/4δ or 1/(δ-δ_{m,n})
    else
        return 1/(P^2 - Pmn^2) # 1/(δ-δ_{m,n})
    end
end

"""
    H(q, Nmax, block, corr, lr;
      der = false, reg = false)

Compute the function ``H(q,δ)``. If der=true, compute instead the function ``H^{\\text{der}}``. If reg=true, compute instead ``H^{\\text{reg}}``.
"""
function H(q, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr;
           der = false, reg = false)
    @assert !(der && reg) "you should not compute the derivative of a regularised block"
    V = block.channelField
    P = V["P"][lr]
    c = corr.charge
    β = c["β"]
    pow = 1

    res = der ? 2*P*log(16*q) : 1 # H_P = 1 + sum(...), H_P^der = 2P log(16q) + sum(...)

    for N in 1:Nmax
        sum_mn = sum(sum(computeCNmn(N, m, n, corr, "s", lr)*block_recursion_coeff(q, c, V, m, n, der, reg, lr)
                         for n in 1:N if m*n <= N) for m in 1:N)

        pow *= 16*q
        res += pow * sum_mn
    end

    return res
end

"""
    block_chiral_schan_value(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

Compute the chiral conformal block

``\\mathcal F^{(s)}_{\\delta}(x)``

"""
function block_chiral_schan(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr;
                            der=false, reg=false)
    return blockprefactor(block, corr, x, lr) * H(qfromx(x), Nmax, block, corr, lr, der=der, reg=reg)
end

"""
    block_chiral(x, Nmax, block, corr, lr)

Compute the chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x)``

where `chan` is `s`, `t`, or `u`."""
function block_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr;
                      der = false, reg = false)
    chan = block.channel
    x_lr = (lr == left ? x : conj(x))
    return channelprefactor_chiral(block, corr, x_lr)*block_chiral_schan(crossratio(chan, x), Nmax, block, permute_ext_fields(corr, chan), lr, der=der, reg=reg)
end

"""
    block_non_chiral(x, Nmax, block, corr)

Compute the non-chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x) \\overline{\\mathcal F}^{(\\text{chan})}_{\\delta}( \bar x )``

where `chan` is `s`,`t` or `u`.

TODO: regularise R_(r,s) / \bar{R}_(r,s)
"""
function block_non_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation)

    chan = block.channel
    Vchan = block.channelField

    if Vchan.isKac && (Vchan.r%1 != 0 || Vchan.s%1 != 0 || spin(Vchan) == 0) # non-logarithmic block

        return channel_sign(block, corr, x) * block_chiral(x, Nmax, block, corr, left) * block_chiral(conj(x), Nmax, block, corr, right)

    elseif 0 == 1 # accidentally non-logarithmic block
        return
    else
        # logarithmic block

        r, s = Vchan.r, Vchan.s

        if Vchan.r < 0 || Vchan.s < 0
            error("Trying to compute a logarithmic block with a negative index: r=$(Vchan.r), s=$(Vchan.s) . This goes against the chosen convention")
        else
            c = corr.charge
            block1 = FourPointBlockSphere(chan, Field(c, Kac=true, r=r, s=s)) # non-log block with momenta (P_(r,s), P_(r,-s)) in the channel
            block2 = FourPointBlockSphere(chan, Field(c, Kac=true, r=r, s=-s)) # non-log block with momenta (P_(r,-s), P_(r,s)) in the channel

            F_Prms = block_chiral(x, Nmax, block2, corr, left) # F_{P_(r,-s)}
            F_Prms_bar = block_chiral(conj(x), Nmax, block1, corr, right) # \bar F_{P_(r,-s)}
            F_der_Prms = block_chiral(x, Nmax, block2, corr, left, der=true) # F'_{P_(r,-s)}
            F_der_Prms_bar = block_chiral(conj(x), Nmax, block1, corr, right, der=true) # \bar F'_{P_(r,-s)}
            F_reg_Prs = block_chiral(x, Nmax, block1, corr, left, reg=true) # F^reg_{P_(r,s)}
            F_reg_Prs_bar = block_chiral(conj(x), Nmax, block2, corr, right, reg=true) # \bar F^reg_{P_(r,s)}

            R = Rmn(r, s, corr, chan, left) # Vchan["P"][left] = P_(r,s)
            R_bar = Rmn(r, s, corr, chan, right)

            term1 = (F_reg_Prs - R*F_der_Prms)*F_Prms_bar
            term2 = R/R_bar*F_Prms*(F_reg_Prs_bar - R_bar*F_der_Prms_bar)
            term3 = -R*ell(corr, r, s)*F_Prms*F_Prms_bar

            return F_Prms, F_Prms_bar, F_der_Prms, F_der_Prms_bar, F_reg_Prs, F_reg_Prs_bar
            # return channel_sign(block, corr, x)*(term1+term2+term3)
        end
    end
end

end # end module

"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

using ..CFTData, ..OnePointCorrelationFunctions
import EllipticFunctions: etaDedekind as η

export OnePointBlockTorus, block

#===========================================================================================
Struct containing the data required to compute a block: an external field
===========================================================================================#
struct OnePointBlockTorus{T}
    channelField::Field{T}
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2

qfromtau(τ) = exp(2im*big(π)*τ)
δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)

#===========================================================================================
Compute the conformal block
===========================================================================================#
"""
    H(q, Nmax, block, corr, leftright)
Compute the function  ``H^{\\text{torus}}(q,δ)``."""
function H(q, Nmax, block::OnePointBlockTorus, corr::OnePointCorrelation, lr)
    δ = block.channelField["δ"][lr]
    B = corr.charge["B"]
    res = 1
    pow = 1
    for N in 1:Nmax
        sum_mn = sum(sum(computeCNmn(N, m, n, corr, lr)/(δ-δrs(m, n, B))
                         for n in 1:N if m*n <= N) for m in 1:N)
        pow *= q
        res += pow * sum_mn
    end
    return res
end

"""
    block_chiral_schan(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

Compute the chiral conformal block

``\\mathcal F^{\text{torus}}_{\\delta}(x)``

"""
function block_chiral(τ, Nmax, block::OnePointBlockTorus, corr::OnePointCorrelation, lr)
    δ = block.channelField["δ"][lr]
    return q^δ/η(τ) * H(qfromtau(τ), Nmax, block, corr, lr)
end

"""
Compute the non-chiral conformal block

`` \\mathcal F_{\\Delta}^{(\\text{chan})}(\\Delta_i| x)``

where ``\\text{chan}`` is `s`,`t` or `u`.

TODO: logarithmic blocks
"""
function F_one_point_torus(τ, Nmax, block::OnePointBlockTorus, corr::OnePointCorrelation)
    block_chiral(τ, Nmax, block, corr, left) * conj(block_chiral(conj(τ), Nmax, block, corr, right))
end

end # end module
