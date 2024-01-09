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
using Match, .CFTData

#export fourpointblock_sphere

#=
Get t- and u- channel blocks from s-channel block
=#
"""Give total prefactor to get t- or u-channel blocks from s-channel block"""
prefactor(block, x) = @match block.channel begin
        "s" => 1
        "t" => (-1)^(sum(spin(block.extFields)))
        "u" => (-1)^(sum(spin.(block.extFields)))*abs2(x)^(-2*block.extFields[1]["Δ"])
    end

"""Give cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
crossratio(channel, x) = @match channel begin
        "s" => x
        "t" => 1-x
        "u" => 1/x
    end

"""Permute the external fields to get t- or u-channel blocks from s-channel block"""
permutationΔext(channel,Δs) = @match channel begin
        "s" => [Δs[1],Δs[2],Δs[3],Δs[4]]
        "t" => [Δs[1],Δs[4],Δs[3],Δs[2]]
        "u" => [Δs[1],Δs[3],Δs[2],Δs[4]]
end;

end

"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

#export onepointblock_torus

end