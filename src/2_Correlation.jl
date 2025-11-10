"""
Abstract supertype for all correlation types.

# Hierarchy

```text
Correlation  (Abstract Type)
├─ Correlation1  (Abstract Type)
│  ├─ ChiralCorrelation1
│  └─ NonChiralCorrelation1
└─ Correlation4  (Abstract Type)
   ├─ ChiralCorrelation4
   └─ NonChiralCorrelation4
```

# Constructors

        Correlation(fields, Δmax)

Create a Correlation object. The fields can be passed as independent arguments or as collections, there can be 1 or 4 of them. The `Correlation`
object stores residues ``R_{m, n}`` and ``C^N_{m, n}`` up to `Δmax`.

Aliased to Corr, Co.

# Examples

```jldoctest; output=false
c = CC(β = 0.5)
V = Field(c, r=2, s=1//2)
co1 = Correlation(V, V, V, V, 10) # NonChiralCorrelation4
co2 = Corr(V, 10) # NonChiralCorrelation1
co3 = Corr(Tuple(V[:left] for _ in 1:4), 20) # ChiralCorrelation4
co3.fields[1] == V[:left]

# output

true
```
"""
abstract type Correlation{T} end
const Corr = Correlation # alias
const Co = Correlation # alias

struct RmnTable{T}
    values::Array{T,2}
    keys::Set{Tuple{Int,Int}}
end

struct CNmnTable{T}
    values::Array{T,3}
    keys::Vector{Set{Tuple{Int,Int}}}
    δs::Matrix{T}
end

abstract type ChiralCorrelation{T} <: Correlation{T} end
const CCo{T} = ChiralCorrelation{T}

struct ChiralCorrelation4{T} <: ChiralCorrelation{T}
    fields::NTuple{4,CD{T}}
    c::CentralCharge{T}
    Rmn::Channels{RmnTable{T}}
    Rmnreg::Channels{RmnTable{T}}
    CNmn::Channels{CNmnTable{T}}
    Δmax::Int
end

struct ChiralCorrelation1{T} <: ChiralCorrelation{T} # N = number of fields
    fields::NTuple{1,CD{T}}
    c::CentralCharge{T}
    Rmn::RmnTable{T}
    Rmnreg::RmnTable{T}
    CNmn::CNmnTable{T}
    Δmax::Int
end

abstract type NonChiralCorrelation{T} <: Correlation{T} end
const NCCo{T} = NonChiralCorrelation{T}

struct NonChiralCorrelation4{T} <: NonChiralCorrelation{T}
    fields::NTuple{4,Field{T}}
    c::CentralCharge{T}
    Rmn::LR{Channels{RmnTable{T}}}
    Rmnreg::LR{Channels{RmnTable{T}}}
    CNmn::LR{Channels{CNmnTable{T}}}
    Δmax::Int
end

struct NonChiralCorrelation1{T} <: NonChiralCorrelation{T}
    fields::NTuple{1,Field{T}}
    c::CentralCharge{T}
    Rmn::LR{RmnTable{T}}
    Rmnreg::LR{RmnTable{T}}
    CNmn::LR{CNmnTable{T}}
    Δmax::Int
end

#===========================================================================================
Four-point residues
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
            s[i] + pm * s[j] isa Real &&
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
            -signs[2] .* (d[i].r, d[i].s) .- signs[2] * signs[1] .* (d[j].r, d[j].s)
        P = ConformalDimension(d[1].c, r = r0, s = s0).P
        return 8 * signs[1] * signs[2] * d[i].P * d[j].P * P
    end
    return 2d[j].P
end

function Rmn_term(r, s, i, j, d::NTuple{4,CD})
    Rmn_term_vanishes(r, s, i, j, d) && return Rmn_term_reg(r, s, i, j, d)
    return Rmn_term_nonzero(r, s, i, j, d)
end

function Rmn_term(r, s, d::NTuple{4,CD{T}})::T where {T}
    prod(Rmn_term(r, s, i, j, d) for (i, j) in ((1, 2), (3, 4)))
end

# write Rmn = \prod_r Pn(r), Pn(r) = \prod_{s=-n+1}^{n-1} Rmn_term(r, s)
function computePns(factors::Matrix{T}, Nmax, _::NTuple{4,CD}) where {T}
    # Pns[n, r+1] = P_n(r), r>=0, n>0
    # factors[(r, s)] = Rmn_term(r, s)
    Pns = Matrix{T}(undef, Nmax, Nmax)
    @inbounds for r = 1:Nmax
        Pns[1, r] = factors[r, 0+Nmax]
        if r - 1 == 0
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
function computeRmns(Nmax, ds::NTuple{4,CD{T}}) where {T}
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
One-point residues
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

Rmn_term_reg(_, _, d::NTuple{1,CD}) = -2 * d[1].P

function Rmn_term(r, s, d::NTuple{1,CD})
    Rmn_term_vanishes(r, s, d) && return Rmn_term_reg(r, s, d)
    return Rmn_term_nonzero(r, s, d)
end

#=
write Rmn = \prod_r Pn(r), Pn(r) = \prod_{s=-n+1}^{n-1} Rmn_term(r, s)
=#
function computePns(factors::Matrix{T}, Nmax, _::NTuple{1,CD}) where {T}
    # Pns[n, r+1] = P_n(r), r>=0, n>0
    Pns = Matrix{T}(undef, Nmax, Nmax)
    for a = 1:Nmax
        r = 2a - 1
        Pns[1, a] = factors[a, Nmax] * factors[a, Nmax+1]
        n = 2
        for n = 2:Nmax
            r * abs(2n - 1) >= 8Nmax && continue
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

function computeRmns(Nmax, ds::NTuple{1,CD{T}}) where {T}
    factors = Matrix{T}(undef, Nmax, 2Nmax)
    for a = 1:Nmax
        for b = 1:2Nmax
            r, s = 2a - 1, 2b - 2Nmax - 1
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
Coefficients CNmn for Zamolodchikov's recursion.
===========================================================================================#
function computeCNmns(Nmax, c::CC{T}, Rs) where {T}
    B = c.B
    mns = [Set{Tuple{Int,Int}}() for _ = 1:Nmax]
    δs = [δrs(mp, np, B) for mp = 1:Nmax, np = 1:Nmax]
    Cs = CNmnTable{T}(Array{T}(undef, (Nmax, Nmax, Nmax)), mns, δs)
    buf = zero(T)
    @inbounds for N = 1:Nmax
        @inbounds for m = 1:Nmax
            maxn = N ÷ m
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
                                # we do
                                # _sum += Cs[Np, mp, np] / (δ0 - Cs.δs[mp, np])
                                # below is the version with mutable arithmetics. much faster because it does not allocate.
                                buf = MA.operate_to!!(buf, -, δ0, Cs.δs[mp, np])
                                buf = MA.operate_to!!(buf, /, Cs[Np, mp, np], buf)
                                _sum = MA.operate_to!!(_sum, +, _sum, buf)
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

"""
   @channels expr(chan)

Create something in each channel.

# Example:
```julia
@channels Block(co, chan, V, 10)
```
expands to
Channels(
    Block(co, :s, V, 10),
    Block(co, :t, V, 10),
    Block(co, :u, V, 10)
)
a second argument can be passed to give a different name
for the chan variable, e.g. this is equivalent to the above
@channels Block(co, myvar, V, 10) myvar"""
macro channels(body, chan = :chan)
    esc(
        quote
            Channels(
                $([:(
                    let $(chan) = $(QuoteNode(name))
                        $(body)
                    end
                ) for name in (:s, :t, :u)]...),
            )
        end,
    )
end

Base.getindex(R::RmnTable, m::Int, n::Int) = R.values[m, n]
function Base.setindex!(R::RmnTable, value, m::Int, n::Int)
    R.values[m, n] = value
    push!(R.keys, (m, n))
    return value
end
Base.haskey(R::RmnTable, key::Tuple{Int,Int}) = key in R.keys

Base.getindex(C::CNmnTable, N::Int, m::Int, n::Int) = C.values[N, m, n]
Base.haskey(C::CNmnTable, k::Tuple{Int,Int,Int}) =
    k[1] < length(C.keys) && (k[2], k[3]) in C.keys[k[1]]
function Base.setindex!(C::CNmnTable, value, N, m::Int, n::Int)
    C.values[N, m, n] = value
    push!(C.keys[N], (m, n))
    return value
end

function permute_4(ds::NTuple{4,T}, chan::Symbol) where {T}
    chan === :s && return ds
    chan === :t && return (ds[1], ds[4], ds[3], ds[2])
    chan === :u && return (ds[1], ds[3], ds[2], ds[4])
end

#======================================================================================
Chiral correlations
======================================================================================#
function ChiralCorrelation(d::NTuple{4,CD{T}}, Δmax::Int) where {T}
    @assert all((dim.c === d[1].c for dim in d)) "
        External fields in the argument of the Correlation constructor do not all have the same
        CentralCharge
    "
    channel_syms = (:s, :t, :u)

    # Launch a task for each channel
    futures = map(1:3) do i
        Threads.@spawn begin
            ch = channel_syms[i]
            dx = permute_4(d, ch)
            r, rreg = computeRmns(Δmax, dx)
            cn = computeCNmns(Δmax, d[1].c, r)
            (r, rreg, cn)
        end
    end
    # Collect results
    vals = Tuple(fetch(f) for f in futures)
    # Extract and regroup the outputs
    Rmn = Channels(Tuple(vals[i][1] for i = 1:3))
    Rmnreg = Channels(Tuple(Tuple(vals[i][2] for i = 1:3)))
    CNmn = Channels(Tuple(Tuple(vals[i][3] for i = 1:3)))

    ChiralCorrelation4{T}(d, d[1].c, Rmn, Rmnreg, CNmn, Δmax)
end

function ChiralCorrelation(d::NTuple{1,CD{T}}, Δmax::Int) where {T}
    # Launch a task for each channel
    Rmn, Rmnreg = computeRmns(Δmax, d)
    CNmn = computeCNmns(Δmax, d[1].c, Rmn)
    ChiralCorrelation1{T}(d, d[1].c, Rmn, Rmnreg, CNmn, Δmax)
end

ChiralCorrelation(d1::CD, d2::CD, d3::CD, d4::CD, Δmax) =
    ChiralCorrelation((d1, d2, d3, d4), Δmax)
ChiralCorrelation(d::CD, Δmax) = ChiralCorrelation(d, Δmax)
ChiralCorrelation(ds::NTuple{4,CD{T}}, c::CC, Rmn, Rmnreg, CNmn, Δmax) where {T} =
    ChiralCorrelation4{T}(ds, c, Rmn, Rmnreg, CNmn, Δmax)
ChiralCorrelation(d::NTuple{1,CD{T}}, c::CC, Rmn, Rmnreg, CNmn, Δmax) where {T} =
    ChiralCorrelation1{T}(d, c, Rmn, Rmnreg, CNmn, Δmax)

function Base.show(io::IO, ::MIME"text/plain", R::RmnTable{T}) where {T}
    print(io, "RmnTable{$T}(")
    for k in sort([k for k in R.keys])
        print(io, "$k => $(R[k...]), ")
    end
    println(io, ")")
end

function Base.show(io::IO, R::RmnTable{T}) where {T}
    println(io, "RmnTable{$T}(")
    for k in sort([k for k in R.keys])
        println(io, "\t$k => $(R[k...]),")
    end
    println(io, ")")
end

function Base.show(io::IO, C::CNmnTable{T}) where {T}
    println(io, "CNmnTable{$T}(")
    for (N, Cs) in enumerate(C.keys)
        print("$N => [")
        for k in sort([k for k in Cs])
            println(io, "\t$k => $(C[N, k[1], k[2]]),")
        end
        println("]")
    end
    println(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", co::CCo{T}) where {T}
    print(io, "ChiralCorrelation{$T} with external dimensions\n$(co.fields)")
end

function Base.show(io::IO, co::CCo)
    print(io, "< ")
    for d in co.fields
        print(io, d, " ")
    end
    print(io, ">")
end

function Base.show(io::IO, ::MIME"text/plain", co::NCCo{T}) where {T}
    println(io, "NonChiralCorrelation{$T} with external fields")
    print(io, "< ")
    for V in co.fields
        print(io, V, " ")
    end
    print(io, ">")
end

function latexstring(co::NCCo)
    res = L"$$\langle "
    for V in co.fields
        res *= strip(latexstring(V), '$')
    end
    res *= L"\rangle$$"
end

function Base.show(io::IO, ::MIME"text/latex", co::NCCo{T}) where {T}
    print(io, latexstring(co))
end

function Base.show(io::IO, co::NCCo)
    print(io, "< ")
    for V in co.fields
        print(io, V, " ")
    end
    print(io, ">")
end

#======================================================================================
Non Chiral correlations
======================================================================================#
NonChiralCorrelation{T}(Vs::NTuple{4,Field}, c, Rmn, Rmnreg, CNmn, Δmax) where {T} =
    NonChiralCorrelation4{T}(Vs, c, Rmn, Rmnreg, CNmn, Δmax)
NonChiralCorrelation{T}(Vs::NTuple{1,Field}, c, Rmn, Rmnreg, CNmn, Δmax) where {T} =
    NonChiralCorrelation1{T}(Vs, c, Rmn, Rmnreg, CNmn, Δmax)

function NonChiralCorrelation(Vs::Tuple{Vararg{Field{T}}}, Δmax) where {T}
    @assert all((V.c === Vs[1].c for V in Vs)) """
    External fields in the argument of the Correlation constructor do not all have the same
    CentralCharge
    """

    dims = LR(Tuple(V[:left] for V in Vs), Tuple(V.dims.right for V in Vs))
    corr_chiral =
        LR(ChiralCorrelation(dims.left, Δmax), ChiralCorrelation(dims.right, Δmax))
    Rmn = LR(corr_chiral.left.Rmn, corr_chiral.right.Rmn)
    Rmnreg = LR(corr_chiral.left.Rmnreg, corr_chiral.right.Rmnreg)
    CNmn = LR(corr_chiral.left.CNmn, corr_chiral.right.CNmn)
    NonChiralCorrelation{T}(Vs, Vs[1].c, Rmn, Rmnreg, CNmn, Δmax)
end

function Base.getindex(c::NCCo, s::Symbol)
    s in (:left, :right) && return ChiralCorrelation(
        Tuple(v.dims[s] for v in c.fields),
        c.c,
        c.Rmn[s],
        c.Rmnreg[s],
        c.CNmn[s],
        c.Δmax,
    )
    error("Only co[:left] and co[:right] are defined")
end

function NonChiralCorrelation(cl::CCo{T}, cr::CCo{T}) where {T}
    dims_left = cl.fields
    dims_right = cr.fields
    fields = Tuple(Field(d1, d2) for (d1, d2) in zip(dims_left, dims_right))
    Δmax = cl.Δmax
    Rmn = LR(cl.Rmn, cr.Rmn)
    Rmnreg = LR(cl.Rmnreg, cr.Rmnreg)
    CNmn = LR(cl.CNmn, cr.CNmn)
    NonChiralCorrelation{T}(fields, cl.c, Rmn, Rmnreg, CNmn, Δmax)
end

function Correlation()
    NonChiralCorrelation(Field(), 0)
end

Correlation(ds::Tuple{Vararg{CD}}, Δmax::Int) = CCo(ds, Δmax)
Correlation(ds::Vector{CD{T}}, Δmax::Int) where {T} = CCo(Tuple(ds), Δmax)
Correlation(Vs::Tuple{Vararg{Field}}, Δmax::Int) = NCCo(Vs, Δmax)
Correlation(Vs::Vector{Field{T}}, Δmax::Int) where {T} = NCCo(Tuple(Vs), Δmax)
Correlation(d1::CD, d2, d3, d4, Δmax::Int) = CCo((d1, d2, d3, d4), Δmax)
Correlation(V1::Field, V2, V3, V4, Δmax::Int) = NCCo((V1, V2, V3, V4), Δmax)
Correlation(d::ConformalDimension, Δmax::Int) = CCo((d,), Δmax)
Correlation(V::Field, Δmax::Int) = NCCo((V,), Δmax)
Correlation(cl::CCo, cr::CCo) = NCCo(cl, cr)

const Correlation4{T} = Union{ChiralCorrelation4{T},NonChiralCorrelation4{T}}
const Correlation1{T} = Union{ChiralCorrelation1{T},NonChiralCorrelation1{T}}
getRmn(co::Correlation4, chan::Symbol) = getfield(co.Rmn, chan)
getRmn(co::Correlation1, _::Symbol) = co.Rmn
getRmnreg(co::Correlation4, chan::Symbol) = getfield(co.Rmnreg, chan)
getRmnreg(co::Correlation1, _::Symbol) = co.Rmnreg
getCNmn(co::Correlation4, chan::Symbol) = getfield(co.CNmn, chan)
getCNmn(co::Correlation1, _::Symbol) = co.CNmn
getfields(co::Correlation4, chan::Symbol) = permute_4(co.fields, chan)
getfields(co::Correlation1, _::Symbol) = co.fields
