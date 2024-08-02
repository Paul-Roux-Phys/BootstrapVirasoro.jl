#=
# This file computes the residues of conformal blocks for four-point sphere correlations,
# for use in the Zamolodchikov recursion.
=#

const channels = (:s, :t, :u)

"""
FourPointCorrelation{T}

Struct for holding data common to all conformal blocks appearing in a given correlation
function:
* Four fields
* Residues `R_mn` in all channels
* Coefficients `C^N_{mn}` in all channels
"""
struct FourPointCorrelation{T} <: Correlation{T}

    fields::NTuple{4, Field{T}}
    Nmax::Int
    _Rmn::Dict{Symbol, RmnTable{T}}
    _Rmn_reg::Dict{Symbol, RmnTable{T}}
    _CNmn::Dict{Symbol, CNmnTable{T}}

end

function Base.show(io::IO, corr::FourPointCorrelation)
    print("Four-point correlation function: < V_1 V_2 V_3 V_4 > where ")
    print("\nV_1 = "); show(corr.fields[1])
    print("\nV_2 = "); show(corr.fields[2])
    print("\nV_3 = "); show(corr.fields[3])
    print("\nV_4 = "); show(corr.fields[4])
end

function permute_fields(Vs::NTuple{4, Field}, chan::Symbol)
    chan === :s && return Vs
    chan === :t && return (Vs[1], Vs[4], Vs[3], Vs[2])
    chan === :u && return (Vs[1], Vs[3], Vs[2], Vs[4])
end

"""TODO: bounds to compute Rmn and CNmn
    FourPointCorrelation(
        c::CentralCharge{T},
        Vs::Tuple{Field{T}, Field{T}, Field{T}, Field{T}},
        Nmax::Int
        )

Constructor function for the FourPointCorrelation type.

# Examples

```julia-repl
julia> using JuliVirBootstrap;
julia> c = CentralCharge(:c, 0.5);
julia> V1 = Field(c, Kac=true, r=1, s=0);
julia> V2 = Field(c, Kac=true, r=2, s=0);
julia> corr = FourPointCorrelation(c, (V1, V2, V2, V1), 6);
julia> corr._Rmn[:s]
Dict{Tuple{Int64, Int64}, Tuple{Float64, Float64}} with 10 entries:
  (3, 2) => (-3.11111e-7, -3.11111e-7)
  (1, 2) => (0.0125355, 0.0125355)
  (3, 1) => (-0.00670399, -0.00670399)
  (1, 1) => (0.158203, 0.158203)
  (5, 1) => (-1.17118e-5, -1.17118e-5)
  (1, 3) => (0.00045679, 0.00045679)
  (2, 2) => (0.00114524, 0.00114524)
  (1, 6) => (1.06026e-7, 1.06026e-7)
  (1, 4) => (2.58811e-5, 2.58811e-5)
  (1, 5) => (1.62729e-6, 1.62729e-6)
```

Note that the parentheses are optional around the four fields:
```julia-repl
julia> corr = FourPointCorrelation(c, V1, V1, V2, V1, 6);
julia> corr._Rmn[:t]
Dict{Tuple{Int64, Int64}, Tuple{Float64, Float64}} with 5 entries:
  (3, 2) => (2.6927e-7, 2.6927e-7)
  (1, 2) => (0.00881619, 0.00881619)
  (2, 2) => (3.11683e-5, 3.11683e-5)
  (1, 6) => (1.59168e-7, 1.59168e-7)
  (1, 4) => (3.32838e-5, 3.32838e-5)
```
"""
function FourPointCorrelation(
    c::CentralCharge{T},
    Vs::NTuple{4, Field{T}},
    Nmax::Int
    ) where {T}
    tmp = Dict(
        s => FourPointCorrelation{T}(permute_fields(Vs, s), Nmax,
            Dict{Symbol,RmnTable{T}}(),
            Dict{Symbol,RmnTable{T}}(),
            Dict{Symbol,CNmnTable{T}}()
        )
        for s in channels
    )

    Rmn_dict = Dict(
        x => Dict((m, n) => (Rmn(m, n, c, tmp[x], left), Rmn(m, n, c, tmp[x], right)) 
                 for m in 1:Nmax, n in 1:Nmax 
                 if Rmn_zero_order(m, n, tmp[x]) == 0 && m*n <= Nmax
                )
        for x in channels
    )

    Rmn_reg = Dict{Symbol, RmnTable{T}}()

    tmp = FourPointCorrelation{T}(Vs, Nmax, Rmn_dict, Rmn_reg, Dict{Symbol, CNmnTable{T}}())
 
    CNmn = Dict(
        x => Dict(
            (N, m, n) => (computeCNmn(N, m, n, c, tmp, x, left),
                          computeCNmn(N, m, n, c, tmp, x, right)
                         )
            for m in 1:Nmax, n in 1:Nmax, N in 1:Nmax
            if m * n <= N && (m, n) in keys(Rmn_dict[x])
        )
        for x in channels
    )

    for x in channels
        filter!(p -> p.second != (0, 0), CNmn[x]) # keep only non-zero values
    end

    FourPointCorrelation{T}(Vs, Nmax, Rmn_dict, Rmn_reg, CNmn)
end

function FourPointCorrelation(
    c::CentralCharge,
    V1::Field, V2::Field, V3::Field, V4::Field,
    Nmax::Int)
    FourPointCorrelation(c, (V1, V2, V3, V4), Nmax)
end


"""
Order of a zero of Rmn, assuming the central charge is generic.
Also return the indices of the vanishing term.
"""
function Rmn_zero_order(m, n, corr::FourPointCorrelation)
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
    for pm in (-1, 1), (i, j) in ((1, 2), (3, 4))
        if (V[i].isKac && V[j].isKac
            && (abs(r[i] + pm * r[j]) <= m - 1 && (r[i] + pm * r[j] - (m - 1)) % 2 == 0)
            && (abs(s[i] + pm * s[j]) <= n - 1 && (s[i] + pm * s[j] - (n - 1)) % 2 == 0))
            
            order += 1

        end
    end

    return order
end

"""Compute one of the terms in the double product of Rmn"""
function Rmn_term(r, s, c::CentralCharge, corr::FourPointCorrelation, lr)
    B = c.B
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
function Rmn(m, n, c::CentralCharge{T}, corr::FourPointCorrelation{T}, lr) where {T}
    if Rmn_zero_order(m, n, corr) == 0
        if m == 1
            res = prod(Rmn_term(0, s, c, corr, lr) for s in 1-n:2:0)
        else # m > 1
            res = prod(prod(Rmn_term(r, s, c, corr, lr)
                            for s in 1-n:2:n-1) for r in 1-m:2:-1)
            if m%2 == 1 # m odd -> treat r=0 term separately
                res *= prod(Rmn_term(0, s, c, corr, lr) for s in 1-n:2:0)
            end
        end
    else
        if m == 1
            res = zero(T)
        end
    end

    return res/(2*Dmn(m, n, c.B))
end

function Rmn_reg(m, n, c::CentralCharge, corr::FourPointCorrelation)
    
end

@memoize function computeCNmn(
    N, m, n,
    c::CentralCharge{T},
    corr::FourPointCorrelation{T}, chan::Symbol, lr
) where {T}
    B = c.B
    rmn = corr._Rmn[chan]
    if !((m, n) in keys(rmn))
        return zero(T)
    elseif m*n > N
        return zero(T)
    elseif m*n == N
        return corr._Rmn[chan][(m, n)][lr]
    else
        res = sum(sum(computeCNmn(N-m*n, mp, np, c, corr, chan, lr)/(δrs(m, -n, B) - δrs(mp, np, B))
                      for mp in 1:N-m*n if mp*np <= N-m*n)
                  for np in 1:N-m*n)
        return corr._Rmn[chan][(m, n)][lr] * res
    end
end