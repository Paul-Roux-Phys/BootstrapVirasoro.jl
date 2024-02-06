#===========================================================================================

ConformalBlocks.jl contains a module ConformalBlocks that computes series expansions for 
Virasoro four-point conformal blocks on the sphere and Virasoro one-point conformal blocks 
on the torus

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram 
Nivesvivat

===========================================================================================#



"""
Series expansion of four-point blocks on the sphere.
"""
module FourPointBlocksSphere

export FourPointBlockSphere, F_four_point_sphere

using ..CFTData, ..FourPointCorrelationFunctions
using Match, EllipticFunctions, Memoization


#===========================================================================================
Struct FourPointBlockSphere
===========================================================================================#
"""
    FourPointBlockSphere{T}

Composite type that represents the list of arguments of a four-point conformal block:
a channel, a field propagating in this channel and four external fields.
"""
struct FourPointBlockSphere{T}

    channel::String
    channelField::Field{T}
    extFields::FourPointCorrelation{T}

end

"""Display blocks"""
function Base.show(io::IO, block::FourPointBlockSphere)
    println("Four-point block")
    println("Channel:\t$(block.channel)")
    println("Channel Field:")
    show(block.channelField)
    println("External Fields:")
    print("1. "); show(block.extFields[1])
    print("2. "); show(block.extFields[2])
    print("3. "); show(block.extFields[3])
    print("4. "); show(block.extFields[4])
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2


#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#
"""Prefactor to get t- or u-channel blocks from the s-channel block"""
function channelprefactor(block::FourPointBlockSphere, x)
    @match block.channel begin
        "s" => 1
        "t" => (-1)^(sum(spin(block.extFields)))
        "u" => (-1)^(sum(spin.(block.extFields)))*abs2(x)^(-2*block.extFields[1]["Δ"])
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

"""Permute the external fields to get t- or u-channel blocks from s-channel block"""
function permute_ext_fields(block::FourPointBlockSphere)
    Vs = block.extFields
    Vs = @match block.channel begin
        "s" => [Vs[1],Vs[2],Vs[3],Vs[4]]
        "t" => [Vs[1],Vs[4],Vs[3],Vs[2]]
        "u" => [Vs[1],Vs[3],Vs[2],Vs[4]]
    end

    return FourPointBlockSphere(block.channel, block.channelField, Vs)
end

#===========================================================================================
Set prefactors, relate the cross-ratio x and the elliptic nome q
===========================================================================================#
"""Nome `q` from the cross-ratio `x`"""
qfromx(x) = exp(-π*ellipticK(1-x) / ellipticK(x))

""""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""Prefactor for getting the block F from H. The argument `lr` indicates if we are working
with a left or right moving block"""
function blockprefactor(block::FourPointBlockSphere, charge::CentralCharge, x, lr)

    e0 = - block.extFields[1]["δ"][lr] - block.extFields[2]["δ"][lr] - (charge["c"]-1)/24
    e1 = - block.extFields[1]["δ"][lr] - block.extFields[4]["δ"][lr] - (charge["c"]-1)/24
    e2 = sum(block.extFields[i]["δ"][lr] for i in 1:4) + (charge["c"]-1)/24
    q=qfromx(x)

    return x^e0 * (1-x)^e1 * jtheta3(0,q)^(-4*e2) * (16*q)^block.channelField[δ][1]
end

#===========================================================================================
Compute the conformal block
===========================================================================================#
function H(q, Nmax, block::FourPointBlockSphere, lr)
    res=0
    δ = block.channelField["δ"][lr]
    for N in 1:Nmax
        for m in 1:N
            for n in 1:N
                if m*n <= N && computeCNmn(N, m, n, B, block, block.channel, lr) != 0
                    res += computeCNmn(N, m, n, B, block, lr) * (16*q)^N/(δ-δrs(m, n, B))
                end
            end
        end
    end
end

"""Compute the chiral conformal block 

``\\mathcal F^{(s)}_{\\delta}(x)``

"""
function Fs_chiral(block::FourPointBlockSphere, charge::CentralCharge, x, lr)
    blockprefactor(block, charge, x, lr) * H(qfromx(x), Nmax, block, lr) 
end

"""Compute the chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x)`` 

where `chan` is `s`, `t`, or `u`."""
function F_chiral(block::FourPointBlockSphere, charge::CentralCharge, x, lr)
    Fs_chiral(permute_ext_fields(block), charge, crossratio(block.channel, x), lr)
end

"""
Compute the non-chiral conformal block

``\\mathcal F_{\\Delta}^{(\\text{chan})}(\\Delta_i| x)``

where `chan` is `s`,`t` or `u`.

TODO: logarithmic blocks
"""
function F_four_point_sphere(block::FourPointBlockSphere, charge::CentralCharge, x)
    channelprefactor(block, x) * Fchiral(block, charge, x, left) * \
    conj(Fchiral(block, charge, conj(x), right))
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
    res=0
    for N in 1:Nmax
        for m in 1:N
            for n in 1:N
                if m*n <= N && computeCNmn(N, m, n, B, block.extField, lr)
                end
            end
        end
    end
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

end #end module