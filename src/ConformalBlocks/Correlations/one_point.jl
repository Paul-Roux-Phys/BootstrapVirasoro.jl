#=
# This file computes the residues of conformal blocks for one-point torus correlations,
# for use in the Zamolodchikov recursion.
=#

"""
OnePointCorrelation{T}

Struct for holding data common to all conformal blocks appearing in a given correlation
function:
* A field object
* Residues `R_mn`
* Coefficients `C^N_{mn}` for the Zamolodchikov recursion.
"""
struct OnePointCorrelation{T} <: Correlation{T}

    field::Field{T}
    _Rmn::RmnTable{T}
    _CNmn::CNmnTable{T}

end

function Base.show(io::IO, corr::OnePointCorrelation)
    println("One-point correlation function: < V > where ")
    print("V = "); show(corr.field)
end

"""TODO: bounds to compute Rmn and CNmn"""
function OnePointCorrelation(c::CentralCharge{T}, V::Field{T}, Nmax::Int) where {T}
    tmp = OnePointCorrelation(V, RmnTable{T}(), CNmnTable{T}())

    Rmn_dict = Dict((m, n) => (Rmn(m, n, c, tmp, left), Rmn(m, n, c, tmp, right)) 
                 for m in 1:Nmax, n in 1:Nmax 
                 if Rmn_zero_order(m, n, tmp) == 0 && m*n <= Nmax
                )

    tmp = OnePointCorrelation(V, Rmn_dict, CNmnTable{T}())
 
    CNmn = Dict(
            (N, m, n) => (computeCNmn(N, m, n, c, tmp, left),
                          computeCNmn(N, m, n, c, tmp, right)
                         )
            for m in 1:Nmax, n in 1:Nmax, N in 1:Nmax
            if m * n <= N && (m, n) in keys(Rmn_dict)
        )

    filter!(p -> p.second != (0, 0), CNmn) # keep only non-zero values

    OnePointCorrelation{T}(V, Rmn_dict, CNmn)
end

"""Order of a pole of Rmn^torus, assuming the central charge is generic"""
function Rmn_zero_order(m, n, corr::OnePointCorrelation)
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
function Rmn(m::Int, n::Int, c::CentralCharge, corr::OnePointCorrelation, lr::Int)
    B = c.B
    V = corr.field
    δ1 = V.δ[lr]
    if Rmn_zero_order(m, n, corr) > 0
        return 0
    else
        res = prod(prod(δrs(r, s, B) - δ1 for r in 1:2:2*m-1) for s in 1-2n:2:2n-1)
        return res/(2*Dmn(m, n, B))
    end
end

@memoize function computeCNmn(
    N::Int,
    m::Int,
    n::Int,
    c::CentralCharge,
    corr::OnePointCorrelation,
    lr::Int
)
    if !((m, n) in keys(corr._Rmn))
        return 0
    elseif m*n > N
        return 0
    elseif m*n == N
        return corr._Rmn[(m, n)][lr]
    else
        B = c.B
        res = sum(sum(computeCNmn(N-m*n, mp, np, c, corr, lr)/(δrs(m, -n, B)-δrs(mp, np, B))
                      for mp in 1:N-m*n if mp*np <= N-m*n)
                  for np in 1:N-m*n)
        return corr._Rmn[(m, n)][lr] * ((N-m*n==0)+res)
    end
end