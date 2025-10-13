include("residues.jl")

struct RmnTable{T}
    values::Array{T,2}
    keys::Set{Tuple{Int,Int}}
end

struct CNmnTable{T}
    values::Array{T,3}
    keys::Vector{Set{Tuple{Int, Int}}}
    δs::Matrix{T}
end

# easily create something in each channel with
# @channels expr(chan)
# example:
# @channels Block(co, chan, V, 10)
# expands to
# Channels(
#     Block(co, :s, V, 10),
#     Block(co, :t, V, 10),
#     Block(co, :u, V, 10)
# )
# a second argument can be passed to give a different name
# for the chan variable, e.g. this is equivalent to the above
# @channels Block(co, myvar, V, 10) myvar
macro channels(body, chan)
    esc(quote
        Channels(
            $( [:( let $(chan) = $(QuoteNode(name))
                        $(body)
                    end ) for name in (:s, :t, :u)]... )
        )
    end)
end
macro channels(body)
    esc(quote
        Channels(
            $( [:( let chan = $(QuoteNode(name))
                        $(body)
                    end ) for name in (:s, :t, :u)]... )
        )
    end)
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

Base.getindex(R::RmnTable, m::Int, n::Int) = R.values[m, n]
function Base.setindex!(R::RmnTable, value, m::Int, n::Int)
    R.values[m, n] = value
    push!(R.keys, (m, n))
    return value
end
Base.haskey(R::RmnTable, key::Tuple{Int,Int}) = key in R.keys

Base.getindex(C::CNmnTable, N::Int, m::Int, n::Int) = C.values[N, m, n]
Base.haskey(C::CNmnTable, k::Tuple{Int,Int,Int}) = k[1] < length(C.keys) && (k[2], k[3]) in C.keys[k[1]]
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

NonChiralCorrelation{T}(Vs::NTuple{4,Field}, c, Rmn, Rmnreg, CNmn, Δmax) where {T} =
    NonChiralCorrelation4{T}(Vs, c, Rmn, Rmnreg, CNmn, Δmax)
NonChiralCorrelation{T}(Vs::NTuple{1,Field}, c, Rmn, Rmnreg, CNmn, Δmax) where {T} =
    NonChiralCorrelation1{T}(Vs, c, Rmn, Rmnreg, CNmn, Δmax)

function NonChiralCorrelation(Vs::Tuple{Vararg{Field{T}}}, Δmax) where {T}
    @assert all((V.c === Vs[1].c for V in Vs)) """
    External fields in the argument of the Correlation constructor do not all have the same
    CentralCharge
    """

    dims = LR(Tuple(V.dims.left for V in Vs), Tuple(V.dims.right for V in Vs))
    corr_chiral = LR(
        ChiralCorrelation(dims.left, Δmax),
        ChiralCorrelation(dims.right, Δmax),
    )
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

const Correlation4{T} = Union{ChiralCorrelation4{T}, NonChiralCorrelation4{T}}
const Correlation1{T} = Union{ChiralCorrelation1{T}, NonChiralCorrelation1{T}}
getRmn(co::Correlation4, chan::Symbol) = getfield(co.Rmn, chan)
getRmn(co::Correlation1, _::Symbol) = co.Rmn
getRmnreg(co::Correlation4, chan::Symbol) = getfield(co.Rmnreg, chan)
getRmnreg(co::Correlation1, _::Symbol) = co.Rmnreg
getCNmn(co::Correlation4, chan::Symbol) = getfield(co.CNmn, chan)
getCNmn(co::Correlation1, _::Symbol) = co.CNmn
getfields(co::Correlation4, chan::Symbol) = permute_4(co.fields, chan)
getfields(co::Correlation1, _::Symbol) = co.fields

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

function Base.show(io::IO, co::NCCo)
    print(io, "< ")
    for V in co.fields
        print(io, V, " ")
    end
    print(io, ">")
end
