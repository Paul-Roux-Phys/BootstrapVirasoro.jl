#===========================================================================================
Four-point case
===========================================================================================#
function double_prod_in_Dmn(m, n, B)
    prod((r * r * B - s * s / B)^2 for s = 1:(n-1) for r = 1:(m-1))
end

function Dmn(m::Int, n::Int, B::T)::T where {T}
    # treat cases m = 1, n=1 separately
    m == 1 && n == 1 && return 1
    m == 1 && return n * prod(s^2 / B * (s^2 / B - m^2 * B) for s = 1:(n-1))
    n == 1 && return m * prod(r^2 * B * (r^2 * B - n^2 / B) for r = 1:(m-1))
    f1 = prod(r^2 * B * (r^2 * B - n^2 / B) for r = 1:(m-1))
    f2 = prod(s^2 / B * (s^2 / B - m^2 * B) for s = 1:(n-1))
    f3 = double_prod_in_Dmn(m, n, B)
    return m * n * f1 * f2 * f3
end

function Rmn_term_vanishes(r, s, i, j, d::NTuple{4,CD})
    !(d[i].isKac && d[j].isKac) && return false

    rs = [d[i].r for i = 1:4]
    ss = [d[i].s for i = 1:4]


    #= The term r, s in Rmn is zero if r1 \pm r2 + r or r3 \pm r4 + r is 0, and
    s1 \pm s2 + s or s3 \pm s4 + s is 0.
    =#
    for pm in (-1, 1), pm2 in (-1, 1)
        if (
            d[i].isKac &&
            d[j].isKac &&
            (rs[i] + pm * rs[j] + pm2 * r == 0) &&
            (ss[i] + pm * ss[j] + pm2 * s == 0)
        )

            return true
        end
    end

    return false
end

function reg_signs(r, s, i, j, d::NTuple{4,CD})
    rs = [d[i].r for i = 1:4]
    ss = [d[i].s for i = 1:4]

    for pm in (-1, 1), pm2 in (-1, 1)
        if (rs[i] + pm * rs[j] + pm2 * r == 0) && (ss[i] + pm * ss[j] + pm2 * s == 0)
            return pm, pm2
        end
    end
end

function Rmn_zero_order(m, n, d::NTuple{4,CD})
    order = 0

    if !((d[1].isKac && d[2].isKac) || (d[3].isKac && d[4].isKac))
        return 0
    end

    r = Tuple(d[i].r for i = 1:4)
    s = Tuple(d[i].s for i = 1:4)

    #= Rmn is zero if r1 \pm r2 or r3 \pm r4 is an integer in 1-m:2:m-1, and
    s1 \pm s2 or s3 \pm s4 is an integer in 1-n:2:n-1.
    equivalently, if (|r1 \pm r2| <= m-1 and r1-r2 - (m-1) % 2 == 0)
    and (|s1 \pm s2| <= n-1 and s1-s2 - (n-1) % 2 == 0)
    =#
    for pm in (-1, 1), (i, j) in ((1, 2), (3, 4))
        if (
            d[i].isKac &&
            d[j].isKac &&
            (
                abs(r[i] + pm * r[j]) <= m - 1 &&
                (r[i] + pm * r[j] - (m - 1)) % 2 == 0
            ) &&
            (abs(s[i] + pm * s[j]) <= n - 1 && (s[i] + pm * s[j] - (n - 1)) % 2 == 0)
        )

            order += 1
        end
    end

    return order
end

function Rmn_term_nonzero(r, s, i, j, d::NTuple{4,CD})
    B = d[1].c.B
    δRS = δrs(r, s, B)
    (r != 0 || s != 0) &&
        return (d[j].δ - d[i].δ)^2 - 2 * δRS * (d[i].δ + d[j].δ) + δRS^2
    return (d[j].δ - d[i].δ) * (-1)^(j / 2)
end

function Rmn_term_reg(r, s, i, j, d::NTuple{4,CD})
    (r != 0 || s != 0) && begin
        signs = reg_signs(r, s, i, j, d)
        r0, s0 =
            -signs[2] .* indices(d[i]) .- signs[2] * signs[1] .* indices(d[j])
        P = ConformalDimension(d[1].c, r = r0, s = s0).P
        return 8 * signs[1] * signs[2] * d[i].P * d[j].P * P
    end
    return 2d[j].P
end

function Rmn_term(r, s, i, j, d::NTuple{4,CD})
    Rmn_term_vanishes(r, s, i, j, d) && return Rmn_term_reg(r, s, i, j, d)
    return Rmn_term_nonzero(r, s, i, j, d)
end

function Rmn_term(r, s, d::NTuple{4, CD{T}})::T where {T}
    prod(Rmn_term(r, s, i, j, d) for (i, j) in ((1, 2), (3, 4)))
end

#=
write Rmn = \prod_r Pn(r), Pn(r) = \prod_{s=-n+1}^{n-1} Rmn_term(r, s)
=#
function computePns(factors::Matrix{T}, Nmax, _::NTuple{4,CD}) where {T}
    # Pns[n, r+1] = P_n(r), r>=0, n>0
    # factors[(r, s)] = Rmn_term(r, s)
    Pns = Matrix{T}(undef, Nmax, Nmax)
    @inbounds for r = 1:Nmax
        Pns[1, r] = factors[r, 0+Nmax]
        if r-1 == 0
            Pns[2, r] = factors[r, 1+Nmax]
        else
            Pns[2, r] = factors[r, 1+Nmax] * factors[r, -1+Nmax]
        end
        @inbounds for n = 3:Nmax
            (r - 1) * (n - 1) > Nmax && break
            Pns[n, r] = Pns[n-2, r] * factors[r, n-1+Nmax]
            if r - 1 != 0
                Pns[n, r] *= factors[r, -n+1+Nmax]
            end
        end
    end
    return Pns
end

function computeDRmns(factors::Matrix{T}, Nmax, ds::NTuple{4,CD}) where {T}
    DRs = Matrix{T}(undef, Nmax, Nmax)
    Pns = computePns(factors, Nmax, ds)
    @inbounds for n = 1:Nmax
        DRs[1, n] = Pns[n, 1]
        DRs[2, n] = Pns[n, 2]
        @inbounds for m = 3:Nmax
            m * n > Nmax && break
            DRs[m, n] = DRs[m-2, n] * Pns[n, m]
        end
    end
    return DRs
end

"""
    computeRmns

Compute the ``Δ``-residues ``R_{m, n}`` for all `m`, `n` such that m*n ≤ Nmax
"""
function computeRmns(Nmax, ds::NTuple{4, CD{T}}) where {T}
    factors = Matrix{T}(undef, Nmax, 2Nmax)
    for r = 1:Nmax
        for s = 1:(2Nmax-1)
            (r - 1) * (s - Nmax) > Nmax && break
            factors[r, s] = Rmn_term(r - 1, s - Nmax, ds)
        end
    end
    DRs = computeDRmns(factors, Nmax, ds)
    Rs = RmnTable{T}(Matrix(undef, Nmax, Nmax), Set{Tuple{Int,Int}}())
    Rregs = RmnTable{T}(Matrix(undef, Nmax, Nmax), Set{Tuple{Int,Int}}())
    for m = 1:Nmax
        for n = 1:Nmax
            m * n > Nmax && break
            if Rmn_zero_order(m, n, ds) == 0
                Rs[m, n] = DRs[m, n] / (2Dmn(m, n, ds[1].c.B))
            else
                Rregs[m, n] = DRs[m, n] / (2Dmn(m, n, ds[1].c.B))
            end
        end
    end
    return Rs, Rregs
end

#===========================================================================================
One-point case
===========================================================================================#
function Rmn_term_vanishes(r, s, D::NTuple{1,CD})
    d = D[1]
    # The term r, s in Rmn is zero if m + r = 0 and n + s = 0
    if d.isKac && d.r + r == 0 && d.s + s == 0
        return true
    end

    return false
end

function Rmn_zero_order(m, n, D::NTuple{1,CD})
    d = D[1]
    if d.isKac &&
       d.r % 2 == 1 &&
       d.s % 2 == 1 &&
       abs(d.r) <= 2 * m - 1 &&
       abs(d.s) <= 2 * n - 1
        return 1
    end
    return 0
end

function Rmn_term_nonzero(r, s, d::NTuple{1,CD})
    B = d[1].c.B
    δ1 = d[1].δ
    return δrs(r, s, B) - δ1
end

Rmn_term_reg(r, s, d::NTuple{1,CD}) = -2 * d[1].P

function Rmn_term(r, s, d::NTuple{1,CD})
    Rmn_term_vanishes(r, s, d) && return Rmn_term_reg(r, s, d)
    return Rmn_term_nonzero(r, s, d)
end

#=
write Rmn = \prod_r Pn(r), Pn(r) = \prod_{s=-n+1}^{n-1} Rmn_term(r, s)
=#
function computePns(factors::Matrix{T}, Nmax, d::NTuple{1,CD}) where {T}
    # Pns[n, r+1] = P_n(r), r>=0, n>0
    Pns = Matrix{T}(undef, Nmax, Nmax)
    for a = 1:Nmax
        r = 2a-1
        Pns[1, a] = factors[a, Nmax] * factors[a, Nmax+1]
        n = 2
        for n = 2:Nmax
            r * abs(2n-1) >= 8Nmax && continue
            Pns[n, a] = Pns[n-1, a] * factors[a, Nmax-n+1] * factors[a, Nmax+n]
        end
    end
    return Pns
end

function computeDRmns(factors::Matrix{T}, Nmax, d::NTuple{1,CD}) where {T}
    Pns = computePns(factors, Nmax, d)
    DRs = Matrix{T}(undef, Nmax, Nmax)
    @inbounds for n = 1:Nmax
        DRs[1, n] = Pns[n, 1]
        @inbounds for m = 2:Nmax
            m * n > Nmax && break
            DRs[m, n] = DRs[m-1, n] * Pns[n, m]
        end
    end
    return DRs
end

function computeRmns(Nmax, ds::NTuple{1, CD{T}}) where {T}
    factors = Matrix{T}(undef, Nmax, 2Nmax)
    for a = 1:Nmax
        for b = 1:2Nmax
            r, s = 2a-1, 2b-2Nmax-1
            r * s >= 8Nmax && continue
            factors[a, b] = Rmn_term(r, s, ds)
        end
    end
    DRs = computeDRmns(factors, Nmax, ds)
    Rs = RmnTable{T}(Matrix(undef, Nmax, Nmax), Set{Tuple{Int,Int}}())
    Rregs = RmnTable{T}(Matrix(undef, Nmax, Nmax), Set{Tuple{Int,Int}}())
    for m = 1:Nmax
        for n = 1:Nmax
            m * n > Nmax && break
            if Rmn_zero_order(m, n, ds) == 0
                Rs[m, n] = DRs[m, n] / (2Dmn(m, n, ds[1].c.B))
            else
                Rregs[m, n] = DRs[m, n] / (2Dmn(m, n, ds[1].c.B))
            end
        end
    end
    return Rs, Rregs
end

#===========================================================================================
Coefficients CNmn
===========================================================================================#
function computeCNmns(Nmax, c::CC{T}, Rs) where {T}
    B = c.B
    mns = [Set{Tuple{Int,Int}}() for _ in 1:Nmax]
    δs = [δrs(mp, np, B) for mp = 1:Nmax, np = 1:Nmax]
    Cs = CNmnTable{T}(Array{T}(undef, (Nmax, Nmax, Nmax)), mns, δs)
    @inbounds for N = 1:Nmax
        @inbounds for m = 1:Nmax
            maxn = N÷m
            @inbounds for n = 1:maxn
                if haskey(Rs, (m, n))
                    Rmn = Rs[m, n]
                    Cs[N, m, n] = Rmn
                    Np = N - m * n
                    Np == 0 && continue
                    δ0 = δrs(m, -n, B)
                    _sum = zero(T)
                    for mp = 1:Np
                        for np = 1:(Np÷mp)
                            if haskey(Cs, (Np, mp, np))
                                _sum += Cs[Np, mp, np] / (δ0 - Cs.δs[mp, np])
                            end
                        end
                    end
                    Cs[N, m, n] *= _sum
                end
            end
        end
    end
    return Cs
end
