#==================

ConformalBlocks.jl contains a module ConformalBlocks that computes series expansions for 
Virasoro four-point conformal blocks on the sphere and Virasoro one-point conformal blocks 
on the torus

==================#


"""
Series expansion of four-point blocks on the sphere
"""
module FourPointBlocksSphere

include("CFTdata.jl")
using Match, .CFTData, EllipticFunctions


#===========================================================================================
Exports
===========================================================================================#
#export fourpointblock_sphere

#===========================================================================================
Struct for the conformal blocks
===========================================================================================#

"""
Object representing a four-point conformal block. The block is represented as the series 
of coefficients of H in Zamolodchikov's recursive formula
"""
struct FourPointBlock{T}

    seriescoefficients::Vector{T}
    channel::String
    prefactor::Function(::T)::T

end

#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#

"""Give total prefactor to get t- or u-channel blocks from s-channel block"""
function channelprefactor(block, x)
    @match block.channel begin
        "s" => 1
        "t" => (-1)^(sum(spin(block.extFields)))
        "u" => (-1)^(sum(spin.(block.extFields)))*abs2(x)^(-2*block.extFields[1]["Δ"])
    end
end

"""Give cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
function crossratio(channel, x)
    @match channel begin
        "s" => x
        "t" => 1-x
        "u" => 1/x
    end
end

"""Permute the external fields to get t- or u-channel blocks from s-channel block"""
function permutationΔext(channel,Δs)
    @match channel begin
        "s" => [Δs[1],Δs[2],Δs[3],Δs[4]]
        "t" => [Δs[1],Δs[4],Δs[3],Δs[2]]
        "u" => [Δs[1],Δs[3],Δs[2],Δs[4]]
    end
end

#===========================================================================================
Set prefactors, relate the various parameters
===========================================================================================#
"""Nome `q` from the cross-ratio `x`"""
function qfromx(x)
    return exp(-π*ellipticK(1-x)/ellipticK(x))
end

"""Prefactor for getting the block F from H"""
function blockprefactor(block, x)
    e1 = block.
end


#===========================================================================================
Implement Zamolodchikov's recursion
===========================================================================================#
end




#=============================================================

=============================================================#




"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

#export onepointblock_torus

end