#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram 
Nivesvivat

===========================================================================================#

module FourPointCorrelationFunctions

using ..CFTData

struct FourPointCorrelation{T}
    Fields::Vector{Field{T}}
end

export FourPointCorrelation, computeCNmn

"""Display a four-point function"""
function Base.show(io::IO, corr::FourPointCorrelation)
    println("Four-point correlation function: < V_1 V_2 V_3 V_4 > where ")
    print("V_1 = "); show(corr.Fields[1])
    print("V_2 = "); show(corr.Fields[2])
    print("V_3 = "); show(corr.Fields[3])
    print("V_4 = "); show(corr.Fields[4])
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2

#===========================================================================================
Compute the coefficients CNmn that are used for computing all conformal blocks in the
expansion of the four-point function
===========================================================================================#
double_prod_in_Dmn(m, n, B) = prod(prod((r^2*B - s^2/B)^2 for s in 1:n-1) for r in 1:m-1)

δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)


function Dmn(m, n, B) 

    if m == 1 && n == 1
        return 1
    elseif m == 1
        return n * prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
    elseif n == 1
        return m * prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
    else
        f1 = prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
        f2 = prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
        f3 = double_prod_in_Dmn(m, n, B)
        return m*n*f1*f2*f3
    end
end

"""Permute the external fields to get t- or u-channel blocks from s-channel block"""
function permute_ext_fields(corr::FourPointCorrelation, channel)
    Vs = corr.Fields
    Vs = @match channel begin
        "s" => [Vs[1],Vs[2],Vs[3],Vs[4]]
        "t" => [Vs[1],Vs[4],Vs[3],Vs[2]]
        "u" => [Vs[1],Vs[3],Vs[2],Vs[4]]
    end

    return FourPointCorrelation(Vs)
end

"""Order of a pole of Rmn, assuming the central charge is generic"""
function Rmn_zero_order(m, n, B, corr::FourPointCorrelation, channel)

    order=0
    V=permute_ext_fields(corr, channel).Fields
    r=[V[i].r for i in 1:4]
    s=[V[i].s for i in 1:4]

    #= Rmn is zero if r1 \pm r2 or r3 \pm r4 is an integer in 1-m:2:m-1, and 
                s1 \pm s2 or s3 \pm s4 is an integer in 1-n:2:n-1.
       equivalently, if (|r1 \pm r2| <= m-1 and r1-r2 - (m-1) % 2 == 0)
                    and (|s1 \pm s2| <= n-1 and s1-s2 - (n-1) % 2 == 0)
    =#
    for pm in (-1,1)
        for (i,j) in ((1,2),(3,4))
            if V[i].isdegenerate && V[j].isdegenerate
                if (abs(r[i]+pm*r[j]) <= m-1 && (r[i]+pm*r[j]-(m-1))%2 == 0) && 
                   (abs(s[i]+pm*s[j]) <= n-1 && (s[i]+pm*s[j]-(n-1))%2 == 0)
                    order += 1
                end
            end
        end
    end

    return order
end

function helper_Rmn(δ1, δ2, δ3, δ4, r, s, B)
    if r == 0 && s == 0
        return (δ2-δ1)*(δ3-δ4)
    else
        return (((δ2-δ1)^2 - 2*δrs(r, s, B)*(δ1+δ2) + δrs(r, s, B)^2)
                *((δ3-δ4)^2 - 2*δrs(r, s, B)*(δ3+δ4) + δrs(r, s, B)^2))
    end
end

"""
Compute `Rmn`. Assume for now there is no singularity.
lr indicates the left or right moving parts of the fields
TODO: value of regularisation
"""
function Rmn(m::Int, n::Int, B, corr::FourPointCorrelation, channel, lr)

    Vs=permute_ext_fields(corr, channel).Fields
    δ1 = Vs[1]["δ"][lr]
    δ2 = Vs[2]["δ"][lr]
    δ3 = Vs[3]["δ"][lr]
    δ4 = Vs[4]["δ"][lr]

    if Rmn_zero_order(m, n, B, corr, channel) > 0
        return 0

    else
        if m == 1
            res = prod(helper_Rmn(δ1, δ2, δ3, δ4, 0, s, B) for s in 1-n:2:0)
        else # m > 1
            res = prod(prod(helper_Rmn(δ1, δ2, δ3, δ4, r, s, B) 
                            for s in 1-n:2:n-1) for r in 1-m:2:-1)
            if m%2 == 1 # m odd -> treat r=0 term separately
                res * prod(helper_Rmn(δ1, δ2, δ3, δ4, 0, s, B) for s in 1-n:2:0)
            end
        end
    end

    return res/(2*Dmn(m, n, B))
end

function computeCNmn(N, m, n, B, corr::FourPointCorrelation, channel, lr)
    res = sum(sum(CNmn(N-m*n, mp, np, B, corr, channel, lr)/(δrs(m, -n, B) - δrs(mp, np, B)) 
                  for mp in 1:N-m*n if mp*np <= N-m*n) 
              for mp in 1:N-m*n)
    return Rmn(m, n, B, corr, channel, lr) * ((N-m*n==0)+res)
end

end # end module