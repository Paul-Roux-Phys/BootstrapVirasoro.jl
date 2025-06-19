const ExtDimensions{T} = Tuple{Vararg{ConformalDimension{T}}} # tuples of dimensions
const FourDimensions{T} = NTuple{4,ConformalDimension{T}}
const OneDimension{T} = Tuple{ConformalDimension{T}}

const ExtFields{T} = Tuple{Vararg{Field{T}}} # tuples of fields
const FourFields{T} = NTuple{4,Field{T}}
const OneField{T} = Tuple{Field{T}}

const ExtPoints{T} = Union{ExtDimensions{T}, ExtFields{T}}
const FourPoints{T} = Union{FourFields{T}, FourDimensions{T}}
const OnePoint{T} = Union{OneField{T}, OneDimension{T}}

include("residues.jl")

"""
    Correlation(args...; Δmax=10.)

Abstract type for holding data relevant to the computation of a correlation function:
* external fields
* Coefficients ``R_{m,n}``, possibly left and right and for different channels
* Coefficients ``R^{\\text{reg}}_{m, n}`` when ``R_{m, n}`` vanishes
* Coefficients ``C^N_{m, n}``, possibly left and right and for different channels
up to `N=Δmax`.

It is also possible to access the central charge, left or right parts of the correlation if it is non chiral, etc. See examples.

# Examples

```jldoctest
julia> c = CentralCharge(β = sqrt(2));

julia> V1 = Field(c, r=2, s=3//2);

julia> co = Correlation(V1, V1, V1, V1, 10)
CorrelationNonChiral{ComplexF64,NTuple{4, Field{ComplexF64}}} with external fields
< V_{(2, 3//2)} V_{(2, 3//2)} V_{(2, 3//2)} V_{(2, 3//2)} >

julia> co._Rmn[:left][:s][2, 4] ≈ 8.37053571428573e-7 - 0.0im
true

julia> co.c.c ≈ -2
true

julia> co.fields
(V_{(2, 3//2)}, V_{(2, 3//2)}, V_{(2, 3//2)}, V_{(2, 3//2)})

julia> co[:left]
CorrelationChiral{ComplexF64} with external dimensions
(Δ_{2, 3//2}, Δ_{2, 3//2}, Δ_{2, 3//2}, Δ_{2, 3//2})
```
"""
abstract type Correlation{T,U} end
# alias
const Corr = Correlation

struct RmnTable{T}
values::Array{T,2}
keys::Set{Tuple{Int,Int}}
end

struct CNmnTable{T}
values::Array{T,3}
keys::Set{Tuple{Int,Int,Int}}
end

const Channels{T} = NamedTuple{(:s, :t, :u),NTuple{3,T}} # type for holding data in all channels

struct CorrelationChiral{T,U<:ExtDimensions{T}} <: Correlation{T,U}
c::CentralCharge{T}
fields::U
Nmax::Int
_Rmn::Channels{RmnTable{T}}
_Rmn_reg::Channels{RmnTable{T}}
_CNmn::Channels{CNmnTable{T}}
end

Base.getindex(R::RmnTable, m::Int, n::Int) = R.values[m, n]
function Base.setindex!(R::RmnTable, value, m::Int, n::Int)
R.values[m, n] = value
push!(R.keys, (m, n))
return value
end
Base.haskey(R::RmnTable, key::Tuple{Int,Int}) = key in R.keys

Base.getindex(C::CNmnTable, N::Int, m::Int, n::Int) = C.values[N, m, n]
Base.haskey(C::CNmnTable, key::Tuple{Int,Int,Int}) = key in C.keys
function Base.setindex!(C::CNmnTable, value, N, m::Int, n::Int)
C.values[N, m, n] = value
push!(C.keys, (N, m, n))
return value
end

function CorrelationChiral(d::U, Nmax::Int) where {T,U<:ExtDimensions{T}}
@assert all((dim.c === d[1].c for dim in d)) """
External fields in the argument of the Correlation constructor do not all have the same
CentralCharge
"""
channel_syms = (:s, :t, :u)
# Launch a task for each channel
futures = map(1:3) do i
    Threads.@spawn begin
        ch = channel_syms[i]
        dx = permute_dimensions(d, ch)
        r, rreg = computeRmns(Nmax, dx)
        cn = computeCNmns!(Nmax, d[1].c, r)
        (r, rreg, cn)
        end
    end

    # Collect results
    vals = Tuple(fetch(f) for f in futures)

    # Now extract and regroup the outputs
    Rmn = Channels{RmnTable{T}}(Tuple(vals[i][1] for i = 1:3))
    Rmnreg = Channels{RmnTable{T}}(Tuple(Tuple(vals[i][2] for i = 1:3)))
    CNmn = Channels{CNmnTable{T}}(Tuple(Tuple(vals[i][3] for i = 1:3)))

    CorrelationChiral{T,U}(d[1].c, d, Nmax, Rmn, Rmnreg, CNmn)
end

CorrelationChiral(d1, d2, d3, d4, Nmax) = CorrelationChiral((d1, d2, d3, d4), Nmax)

struct CorrelationNonChiral{T,U} <: Correlation{T,U}
    c::CentralCharge{T}
    fields::ExtFields{T}
    Nmax::Int
    _Rmn::LeftRight{Channels{RmnTable{T}}}
    _Rmn_reg::LeftRight{Channels{RmnTable{T}}}
    _CNmn::LeftRight{Channels{CNmnTable{T}}}
end

function permute_dimensions(ds::FourDimensions, chan::Symbol)
    chan === :s && return ds
    chan === :t && return (ds[1], ds[4], ds[3], ds[2])
    chan === :u && return (ds[1], ds[3], ds[2], ds[4])
end

permute_dimensions(d::OneDimension, chan::Symbol) = d

function permute_fields(Vs::FourFields, chan::Symbol)
    chan === :s && return Vs
    chan === :t && return (Vs[1], Vs[4], Vs[3], Vs[2])
    chan === :u && return (Vs[1], Vs[3], Vs[2], Vs[4])
end
permute_fields(V::OneField, chan::Symbol) = V

channels(Vs::FourPoints) = (:s, :t, :u)
channels(V::OnePoint) = (:s, :t, :u)

function CorrelationNonChiral(V::U, Nmax::Int) where {T,U<:ExtFields{T}}
    @assert all((v.c === V[1].c for v in V)) """
    External fields in the argument of the Correlation constructor do not all have the same
    CentralCharge
    """

    dims = Tuple(Tuple(v.dims[lr] for v in V) for lr in (:left, :right))
    corr_chiral = Tuple(CorrelationChiral(dims[lr], Nmax) for lr in (:left, :right))

    Rmn = Tuple(corr_chiral[lr]._Rmn for lr in (:left, :right))
    Rmnreg = Tuple(corr_chiral[lr]._Rmn_reg for lr in (:left, :right))
    CNmn = Tuple(corr_chiral[lr]._CNmn for lr in (:left, :right))

    CorrelationNonChiral{T,U}(V[1].c, V, Nmax, Rmn, Rmnreg, CNmn)
end

function Base.getindex(c::CorrelationNonChiral, s::Symbol)
    s in (:left, :right) && return CorrelationChiral(
        c.c,
        Tuple(v.dims[s] for v in c.fields),
        c.Nmax,
        c._Rmn[s],
        c._Rmn_reg[s],
        c._CNmn[s],
    )
    error("cannot access Correlation")
end

function CorrelationNonChiral(corr_chiral::LeftRight{CorrelationChiral{T,U}}) where {T,U}
    dims_left = corr_chiral[:left].fields
    dims_right = corr_chiral[:left].fields
    fields = Tuple(Field(d1, d2) for (d1, d2) in zip(dims_left, dims_right))
    Nmax = corr_chiral[:left].Nmax

    Rmn = Tuple(corr_chiral[lr]._Rmn for lr in (:left, :right))
    Rmnreg = Tuple(corr_chiral[lr]._Rmn_reg for lr in (:left, :right))
    CNmn = Tuple(corr_chiral[lr]._CNmn for lr in (:left, :right))

    CorrelationNonChiral{T,U}(corr_chiral[1].c, fields, Nmax, Rmn, Rmnreg, CNmn)
end

function Correlation()
    CorrelationNonChiral(Field(), 0)
end

Correlation(ds::ExtDimensions, Nmax::Int) = CorrelationChiral(ds, Nmax)
Correlation(Vs::ExtFields, Nmax::Int) = CorrelationNonChiral(Vs, Nmax)
Correlation(Vs::ExtFields, lr::Symbol, Nmax::Int) =
    CorrelationChiral(Tuple(v.dims[lr] for v in Vs), Nmax)
Correlation(d1::ConformalDimension, d2, d3, d4, Nmax::Int) =
    CorrelationChiral((d1, d2, d3, d4), Nmax)
Correlation(V1::Field, V2, V3, V4, Nmax::Int) = CorrelationNonChiral((V1, V2, V3, V4), Nmax)
Correlation(V1::Field, V2, V3, V4, lr::Symbol, Nmax::Int) =
    CorrelationChiral((V1.dims[lr], V2.dims[lr], V3.dims[lr], V4.dims[lr]), Nmax)
Correlation(d::ConformalDimension, Nmax::Int) = CorrelationChiral((d,), Nmax)
Correlation(V::Field, Nmax::Int) = CorrelationNonChiral((V,), Nmax)
Correlation(V::Field, lr, Nmax::Int) = CorrelationChiral((V.dims[lr],), Nmax)
Correlation(co_left::CorrelationChiral, co_right::CorrelationChiral) =
    CorrelationNonChiral((co_left, co_right))
Correlation(args...; Δmax = 10.0) = Correlation(args..., N_max(CD(), Δmax))

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
    for k in sort([k for k in C.keys])
        println(io, "\t$k => $(C[k...]),")
    end
    println(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", co::CorrelationChiral{T}) where {T}
    print(io, "CorrelationChiral{$T} with external dimensions\n$(co.fields)")
end

function Base.show(io::IO, co::CorrelationChiral)
    print(io, "< ")
    for d in co.fields
        print(io, d, " ")
    end
    print(io, ">")
end

function Base.show(io::IO, ::MIME"text/plain", co::CorrelationNonChiral{T,U}) where {T,U}
    println(io, "CorrelationNonChiral{$T,$U} with external fields")
    print(io, "< ")
    for V in co.fields
        print(io, V, " ")
    end
    print(io, ">")
end

function Base.show(io::IO, co::CorrelationNonChiral)
    print(io, "< ")
    for V in co.fields
        print(io, V, " ")
    end
    print(io, ">")
end
