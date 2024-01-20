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

include("CFTdata.jl")

using Match, .CFTData, EllipticFunctions, Memoize


#===========================================================================================
Exports
===========================================================================================#
export FourPointBlockSphere, F_four_point_sphere



#===========================================================================================
Struct containing the data required to compute a block: a channel, four external fields and
a channel field.
===========================================================================================#

struct FourPointBlockSphere{T}

    channel::String
    channelField::Field{T}
    extFields::Vector{Field{T}}

end

function FourPointBlockSphere(channel, channelField::Field{T}, extFields::Vector{Field{T}}) where {T}
    return FourPointBlockSphere{T}(channel, channelField, extFields)
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2 



#===========================================================================================
Compute the conformal block
===========================================================================================#

Fs_chiral(block::FourPointBlockSphere, charge::CentralCharge, x) = \
    blockprefactor(block, charge, x) * H(q, δ, Nmax, block) 

F_chiral(block::FourPointBlockSphere, charge::CentralCharge, x, lr) = \
    Fs_chiral(permute_ext_fields(block), charge, crossratio(block.channel, x))


"""
Compute the non-chiral conformal block
```math
\\mathcal F_{\\Delta}^{(c)}(\\Delta_i| x) where `c` is `s`,`t` or `u`
```
TODO: logarithmic blocks
"""
F_four_point_sphere(block::FourPointBlockSphere, charge::CentralCharge, x) = \
    channelprefactor(block, x) * Fchiral(block, charge, x, left) * \
    conj(Fchiral(block, charge, conj(x), right))



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
function permute_ext_fields(block)
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
Implement Zamolodchikov's recursion
===========================================================================================#
double_prod_in_Dmn(m, n, B) = prod(prod((r^2*B - s^2/B)^2 for s in 1:n-1) for r in 1:m-1)

δrs(r, s, B) = -1/4(B*r^2 + 2r*s + s^2/B)

function Dmn(m, n, B) 

    f1 = prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
    f2 = prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
    f3 = double_prod_in_Dmn(m, n, B)
    
    return m*n*f1*f2*f3
end

"""Order of a pole of Rmn, assuming the central charge is generic"""
function Rmn_pole_order(m, n, B, block::FourPointBlockSphere)

    order=0
    V=block.extFields
    r=[V[i].r for i in 1:4]
    s=[V[i].s for i in 1:4]

    #= pole if r1 \pm r2 or r3 \pm r4 is an integer between 1-m and m-1, and 
                s1 \pm s2 or s3 \pm s4 is an integer between 1-n and n-1.
    =#
    for pm in (-1,1)
        for (i,j) in ((1,2),(3,4))
            if V[i].degenerate && V[j].degenerate
                if r[i]+pm*r[j] in 1-m:2:m-1 && s[i]+pm*s[j] in 1-n:2:n-1
                    order += 1
                end
            end
        end
    end

    return order
end

"""
Compute `Rmn`. Assume for now there is no singularity.
lr indicates a left or right moving block
TODO: value of regularisation
"""
function Rmn(m::Int, n::Int, B, block::FourPointBlockSphere, lr)

    Vs=block.extFields
    δ1 = Vs[1]["δ"][lr]
    δ2 = Vs[2]["δ"][lr]
    δ3 = Vs[3]["δ"][lr]
    δ4 = Vs[4]["δ"][lr]

    if Rmn_pole_order(m, n, B, block) > 0
        return 0

    else
        return ( 1/(2*Dmn) * prod(prod(
                                ((δ2-δ1)^2 - 2δrs(r, s, B)*(δ1+δ2) + δrs(r, s, B)^2)
                               *((δ3-δ4)^2 - 2δrs(r, s, B)*(δ3+δ4) + δrs(r, s, B)^2)
                  for s in 1-n:2:n-1) for r in 1-m:2:m-1)
        )
    end

end

@memoize function CNmn(N, m, n, B, block::FourPointBlockSphere, lr)
    res = sum(sum(CNmn(N-m*n, mp, np, B, block, lr)/(δrs(m, -n, B) - δrs(mp, np, B)) 
                                                        for mp in 1:N-m*n if mp*np <= N-m*n) 
              for mp in 1:N-m*n)
    return Rmn(m, n, B, block, lr) * ((N-m*n == 0) + res)
end

H(q, δ, Nmax, block::FourPointBlockSphere, lr) = \
    1 + sum(sum(sum(CNmn(N, m, n, B, block, lr) * (16q)^N/(δ - δrs(m, n, B))
                    for n in 1:N if m*n <= N) for m in 1:N) for N in 1:Nmax)

end # end module



"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

export F_one_point_torus

#===========================================================================================
Struct containing the data required to compute a block: an external field
===========================================================================================#

struct OnePointBlockTorus{T}

    extField::Field{T}

end

function OnePointBlockTorus(extField::Field{T}) where {T}
    return OnePointBlockTorus{T}(extField)
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2



#===========================================================================================
Compute the conformal block
===========================================================================================#

F_chiral(block::FourPointBlockSphere, charge::CentralCharge, x, lr) = \
    Fs_chiral(permute_ext_fields(block), charge, crossratio(block.channel, x))


"""
Compute the non-chiral conformal block
```math
\\mathcal F_{\\Delta}^{(c)}(\\Delta_i| x) where `c` is `s`,`t` or `u`
```
TODO: logarithmic blocks
"""
F_four_point_sphere(block::FourPointBlockSphere, charge::CentralCharge, x) = \
    channelprefactor(block, x) * Fchiral(block, charge, x, left) * \
    conj(Fchiral(block, charge, conj(x), right))


end #end module