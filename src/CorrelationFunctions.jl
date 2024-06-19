#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module FourPointCorrelationFunctions

export FourPointCorrelation, computeCNmn

using ..CFTData
using Match
import Memoization: @memoize

"""Struct representing a four-point function. Contains
- a central charge
- 4 external fields
"""
struct FourPointCorrelation{T}
    charge::CentralCharge{T}
    fields::Vector{Field{T}}
end

function FourPointCorrelation(c::CentralCharge, V1, V2, V3, V4)
    return FourPointCorrelation(c, [V1, V2, V3, V4])
end

"""Display a four-point function"""
function Base.show(io::IO, corr::FourPointCorrelation)
    print("Four-point correlation function: < V_1 V_2 V_3 V_4 > where ")
    print("\nV_1 = "); show(corr.fields[1])
    print("\nV_2 = "); show(corr.fields[2])
    print("\nV_3 = "); show(corr.fields[3])
    print("\nV_4 = "); show(corr.fields[4])
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2

double_prod_in_Dmn(m, n, B) = prod(prod((r^2*B - s^2/B)^2 for s in 1:n-1) for r in 1:m-1)

δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)

function Dmn(m, n, B)
    if m == 1 && n == 1 # treat cases m = 1, n=1 separately
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

"""Order of a zero of Rmn, assuming the central charge is generic. Also return the indices of the vanishing term."""
function Rmn_zero_order(m, n, corr::FourPointCorrelation)
    B = corr.charge.B
    order = 0
    V=corr.fields

    if !((V[1].isKac && V[2].isKac) || (V[3].isKac && V[4].isKac))
        return 0
    end

    r=[V[i].r for i in 1:4]
    s=[V[i].s for i in 1:4]

    #= Rmn is zero if r1 \pm r2 or r3 \pm r4 is an integer in 1-m:2:m-1, and
    s1 \pm s2 or s3 \pm s4 is an integer in 1-n:2:n-1.
    equivalently, if (|r1 \pm r2| <= m-1 and r1-r2 - (m-1) % 2 == 0)
    and (|s1 \pm s2| <= n-1 and s1-s2 - (n-1) % 2 == 0)
    =#
    for pm in (-1,1)
        for (i,j) in ((1,2), (3,4))
            if V[i].isKac && V[j].isKac
                if (abs(r[i]+pm*r[j]) <= m-1 && (r[i]+pm*r[j]-(m-1))%2 == 0) &&
                    (abs(s[i]+pm*s[j]) <= n-1 && (s[i]+pm*s[j]-(n-1))%2 == 0)
                    order += 1
                end
            end
        end
    end

    return order
end

"""Compute one of the terms in the double product of Rmn"""
function Rmn_term(r, s, corr::FourPointCorrelation, lr)
    B = corr.charge.B
    V = corr.fields
    δ = [V[i].δ[lr] for i in 1:4]
    if r != 0 || s != 0
        return (((δ[2]-δ[1])^2 - 2*δrs(r, s, B)*(δ[1]+δ[2]) + δrs(r, s, B)^2)
                *((δ[3]-δ[4])^2 - 2*δrs(r, s, B)*(δ[3]+δ[4]) + δrs(r, s, B)^2))
    else
        return (δ[2]-δ[1])*(δ[3]-δ[4])
    end
end

"""Compute the regularization of a term in the double product of Rmn"""
function Rmn_term_reg(r, s, corr::FourPointCorrelation, lr)
    V = corr.fields
    if r != 0 || s != 0
        return 8*V[1].P[lr]*V[2].P[lr]*Field(corr.charge, Kac=true, r=r, s=s)
    else
        return 2*V[2].P[lr]
    end
end

"""
Compute `Rmn`.
lr indicates the left or right moving parts of the fields
Cache the result.
TODO: value of regularisation
"""
@memoize function Rmn(m, n, corr::FourPointCorrelation, lr)

    if Rmn_zero_order(m, n, corr) == 0
        if m == 1
            res = prod(Rmn_term(0, s, corr, lr) for s in 1-n:2:0)
        else # m > 1
            res = prod(prod(Rmn_term(r, s, corr, lr)
                            for s in 1-n:2:n-1) for r in 1-m:2:-1)
            if m%2 == 1 # m odd -> treat r=0 term separately
                res *= prod(Rmn_term(0, s, corr, lr) for s in 1-n:2:0)
            end
        end
    else
        if m == 1
            res = 0
        end
    end

    return res/(2*Dmn(m, n, corr.charge.B))
end

@memoize function computeCNmn(N, m, n, corr::FourPointCorrelation, lr)
    B = corr.charge.B
    if Rmn_zero_order(m, n, corr) > 0
        return 0
    elseif m*n > N
        return 0
    elseif m*n == N
        return Rmn(m, n, corr, lr)
    else
        res = sum(sum(computeCNmn(N-m*n, mp, np, corr, lr)/(δrs(m, -n, B) - δrs(mp, np, B))
                      for mp in 1:N-m*n if mp*np <= N-m*n)
                  for np in 1:N-m*n)
        return Rmn(m, n, corr, lr) * res
    end
end

end # end module

module OnePointCorrelationFunctions

export OnePointCorrelation, computeCNmn

using ..CFTData
import ..FourPointCorrelationFunctions: Dmn, δrs # re-use the Dmn from four-point functions

struct OnePointCorrelation{T}
    charge::CentralCharge{T}
    field::Field{T}
end

"""Display a one-point function"""
function Base.show(io::IO, corr::OnePointCorrelation)
    println("One-point correlation function: < V > where ")
    print("V = "); show(corr.field)
end

"""Order of a pole of Rmn^torus, assuming the central charge is generic"""
function Rmn_zero_order(m, n, corr::OnePointCorrelation)
    B = corr.charge.B
    V = corr.field
    if V.isKac && V.r%2==1 && V.s%2==1 && abs(V.r) <= 2*m-1 && abs(V.s) <= 2*n-1
        return 1
    end
    return 0
end

"""
Compute `Rmn^torus`.
lr indicates the left or right moving parts of the fields
TODO: value of regularisation
"""
function Rmn(m, n, corr::OnePointCorrelation, lr)
    B = corr.charge.B
    V = corr.field
    δ1 = V.δ[lr]
    if Rmn_zero_order(m, n, corr) > 0
        return 0
    else
        res = prod(prod(δrs(r, s, B) - δ1 for r in 1:2:2*m-1) for s in 1-2n:2:2n-1)
        return res/(2*Dmn(m, n, B))
    end
end

function computeCNmn(N, m, n, corr::OnePointCorrelation, lr)
    B = corr.charge.B
    if Rmn_zero_order(m, n, corr) > 0
        return 0
    elseif m*n > N
        return 0
    elseif m*n == N
        return Rmn(m, n, corr, lr)
    else
        res = sum(sum(computeCNmn(N-m*n, mp, np, corr, lr)/(δrs(m, -n, B)-δrs(mp, np, B))
                      for mp in 1:N-m*n if mp*np <= N-m*n)
                  for np in 1:N-m*n)
        return Rmn(m, n, corr, lr) * ((N-m*n==0)+res)
    end
end

end # end module
