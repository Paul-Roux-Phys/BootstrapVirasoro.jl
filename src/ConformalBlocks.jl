#===========================================================================================

ConformalBlocks.jl contains a module ConformalBlocks that computes series expansions for 
Virasoro four-point conformal blocks on the sphere and Virasoro one-point conformal blocks 
on the torus

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram 
Nivesvivat

===========================================================================================#


"""
Series expansion of four-point blocks on the sphere
"""
module FourPointBlocksSphere

include("CFTdata.jl")

using Match, .CFTData, EllipticFunctions


#===========================================================================================
Exports
===========================================================================================#
export FourPointBlockSphere

#===========================================================================================
Struct for the conformal blocks
===========================================================================================#

"""
Object representing a four-point conformal block. The block is represented as the series 
of coefficients of H in Zamolodchikov's recursive formula
"""
struct FourPointBlockSphere{T}

    seriescoefficients::Vector{T}
    channel::String
    channelField::Field{T}
    extFields::Vector{Field{T}}
    prefactor::T

end

#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#

"""Prefactor to get t- or u-channel blocks from s-channel block"""
function channelprefactor(block, x)
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
function permutationΔext(channel, Vs)
    @match channel begin
        "s" => [Δs[1],Δs[2],Δs[3],Δs[4]]
        "t" => [Δs[1],Δs[4],Δs[3],Δs[2]]
        "u" => [Δs[1],Δs[3],Δs[2],Δs[4]]
    end
end

#===========================================================================================
Set prefactors, relate the x and the elliptic nome q
===========================================================================================#
"""Nome `q` from the cross-ratio `x`"""
qfromx(x) = exp(-π*ellipticK(1-x) / ellipticK(x))

xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""Prefactor for getting the block F from H"""
function blockprefactor(block, charge, x)

    e0 = - block.extFields[1]["δ"][1] - block.extFields[2]["δ"][1] - (charge["c"]-1)/24
    e1 = - block.extFields[1]["δ"][1] - block.extFields[4]["δ"][1] - (charge["c"]-1)/24
    e2 = sum(block.extFields[i]["δ"][1] for i in 1:4) + (charge["c"]-1)/24

    q=qfromx(x)

    return x^e0 * (1-x)^e1 * jtheta3(0,q)^(-4*e2) * (16*q)^block.channelField[δ][1]
end

#===========================================================================================
Implement Zamolodchikov's recursion
===========================================================================================#
function Rmn(block, )

    for r in n-1:-2:0

    end
end

#===========================================================================================
Block constructor
===========================================================================================#

end




"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

#export onepointblock_torus

end