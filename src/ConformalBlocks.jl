#===========================================================================================

ConformalBlocks.jl contains modules that compute series expansions for
Virasoro four-point conformal blocks on the sphere and Virasoro one-point conformal blocks
on the torus.

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#


"""
Series expansion of four-point blocks on the sphere.
"""
module FourPointBlocksSphere

export FourPointBlockSphere, F_four_point_sphere

using ..CFTData, ..FourPointCorrelationFunctions
import ..FourPointCorrelationFunctions: permute_ext_fields
using Match, EllipticFunctions, Memoization

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
function channelprefactor(block::FourPointBlockSphere, corr::FourPointCorrelation, x)
    @match block.channel begin
        "s" => 1
        "t" => (-1)^(sum(spin(corr.fields)))
        "u" => (-1)^(sum(spin.(corr.fields)))*abs2(x)^(-2*corr.fields[1]["Δ"])
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
@memoize qfromx(x) = exp(-π*ellipticK(1-x) / ellipticK(x))

""""Cross ratio `x` from the nome `q`"""
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

#===========================================================================================
Compute the conformal block
===========================================================================================#
"""Compute the function ``H(q,δ)``."""
function H(q, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr)
    δ = block.channelField["δ"][lr]
    B = corr.charge["B"]
    sq = 16*q
    res=1
    pow = 1
    for N in 1:Nmax
        sum_mn = sum(sum(computeCNmn(N, m, n, corr, block.channel, lr)/(δ-δrs(m, n, B))
                         for n in 1:N if m*n <= N) for m in 1:N)
        # for m in 1:N
        #     for n in 1:N
        #         C = @memoize computeCNmn(N, m, n, corr, block.channel, lr)
        #         if m*n <= N && C != 0
        #             sum_mn += C*pow/(δ-δrs(m, n, B))
        #         end
        #     end
        # end
        pow *= 16*q
        res += pow * sum_mn
    end
    return res
end

"""
    Fs_chiral(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

    Compute the chiral conformal block

``\\mathcal F^{(s)}_{\\delta}(x)``

"""
function Fs_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr)
    blockprefactor(block, corr, x, lr) * H(qfromx(x), Nmax, block, corr, lr)
end

"""Compute the chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x)``

where `chan` is `s`, `t`, or `u`."""
function F_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr)
    chan = block.channel
    Fs_chiral(crossratio(chan, x), Nmax, block, permute_ext_fields(corr, chan), lr)
end

"""
Compute the non-chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x) \\overline{\\mathcal F}^{(\\text{chan})}_{\\delta}( \bar x )``

where `chan` is `s`,`t` or `u`.

TODO: logarithmic blocks
"""
function F_four_point_sphere(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation)
    channelprefactor(block, corr, x) * \
        F_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, left) * \
        conj(F_chiral(conj(x), Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, right))
end

end # end module

"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

using ..CFTData, ..OnePointCorrelationFunctions

export F_one_point_torus

#===========================================================================================
Struct containing the data required to compute a block: an external field
===========================================================================================#
struct OnePointBlockTorus{T}
    extField::OnePointCorrelation{T}
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2

#===========================================================================================
Compute the conformal block
===========================================================================================#
function H(q, Nmax, block::OnePointBlockTorus, lr)
    δ = block.extField["δ"][lr]
    res=1
    prev=0
    next=0
    for N in 1:Nmax
        for m in 1:N
            for n in 1:N
                C = computeCNmn(N, m, n, block.extField, lr)
                if m*n <= N && C != 0
                    res += C*1
                end
            end
        end
    end
    return res
end

function F_chiral(block::OnePointBlockTorus, charge::CentralCharge, x, lr)
end

"""
Compute the non-chiral conformal block

`` \\mathcal F_{\\Delta}^{(\\text{chan})}(\\Delta_i| x)``

where ``\\text{chan}`` is `s`,`t` or `u`.

TODO: logarithmic blocks
"""
function F_one_point_torus(block::OnePointBlockTorus, charge::CentralCharge, x)
    channelprefactor(block, x) * Fchiral(block, charge, x, left) * \
        conj(Fchiral(block, charge, conj(x), right))
end

end # end module
