struct RmnTable{T}
    values::Array{T, 2}
    keys::Set{Tuple{Int, Int}}
end

struct CNmnTable{T}
    values::Array{T, 3}
    keys::Set{Tuple{Int, Int, Int}}
end

const Channels{T} = NamedTuple{(:s, :t, :u), NTuple{3, T}} # type for holding data in all channels

struct CorrelationChiral{T} <: Correlation{T}

    dims::ExtDimensions{T}
    Nmax::Int
    _Rmn::Channels{RmnTable{T}}
    _Rmn_reg::Channels{RmnTable{T}}
    _CNmn::Channels{CNmnTable{T}}

end

function Base.getproperty(co::CorrelationChiral, s::Symbol)
    s === :fields && return co.dims
    getfield(co, s)
end

Base.getindex(R::RmnTable, m::Int, n::Int) = R.values[m, n]
function Base.setindex!(R::RmnTable, value, m::Int, n::Int)
    R.values[m, n] = value
    push!(R.keys, (m, n))
    return value
end
Base.haskey(R::RmnTable, key::Tuple{Int, Int}) = key in R.keys

Base.getindex(C::CNmnTable, N::Int, m::Int, n::Int) = C.values[N, m, n]
Base.haskey(C::CNmnTable, key::Tuple{Int, Int, Int}) = key in C.keys
function Base.setindex!(C::CNmnTable, value, N, m::Int, n::Int)
    C.values[N, m, n] = value
    push!(C.keys, (N, m, n))
    return value
end

function Base.show(io::IO, ::MIME"text/plain", co::CorrelationChiral{T}) where {T}
    print(io, "CorrelationChiral{$T} with external dimensions\n$(co.dims)")
end

function Base.show(io::IO, co::CorrelationChiral)
    print(io, "< ")
    for d in co.dims
        print(io, d, " ")
    end
    print(io, ">")
end

function CorrelationChiral(d::ExtDimensions{T}, Nmax::Int) where {T}
    @assert all((dim.c === d[1].c for dim in d)) """
    External fields in the argument of the Correlation constructor do not all have the same
    CentralCharge
    """
    DRs = Matrix{T}(undef, (Nmax, Nmax))
    Pns = Matrix{T}(undef, (Nmax, Nmax))
    factors = Matrix{T}(undef, (Nmax, 2Nmax))

    channel_syms = (:s, :t, :u)

    vals = Tuple(
        begin
            ch = channel_syms[i]
            dx = permute_dimensions(d, ch)
            r, rreg = computeRmns!(DRs, Pns, factors, Nmax, dx)
            cn = computeCNmns!(Nmax, d[1].c, r)
            r, rreg, cn
        end
        for i in 1:3
    )

    # Now extract and regroup the outputs
    Rmn = Channels{RmnTable{T}}(Tuple(vals[i][1] for i in 1:3))
    Rmnreg = Channels{RmnTable{T}}(Tuple(Tuple(vals[i][2] for i in 1:3)))
    CNmn = Channels{CNmnTable{T}}(Tuple(Tuple(vals[i][3] for i in 1:3)))

    CorrelationChiral{T}(d, Nmax, Rmn, Rmnreg, CNmn)
end

CorrelationChiral(d1, d2, d3, d4, Nmax) = CorrelationChiral((d1, d2, d3, d4), Nmax)

struct CorrelationNonChiral{T} <: Correlation{T}

    fields::ExtFields{T}
    Nmax::Int
    _Rmn::LeftRight{Channels{RmnTable{T}}}
    _Rmn_reg::LeftRight{Channels{RmnTable{T}}}
    _CNmn::LeftRight{Channels{CNmnTable{T}}}

end

function Base.show(io::IO, ::MIME"text/plain", co::CorrelationNonChiral{T}) where {T}
    println(io, "CorrelationNonChiral{$T} with external fields")
    print(io, co.fields)
end

function Base.show(io::IO, co::CorrelationNonChiral)
    print(io, "< ")
    for V in co.fields
        print(io, V, " ")
    end
    print(io, ">")
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

permute_dimensions(V::OneField, chan::Symbol) = V

channels(Vs::FourDimensions) = (:s, :t, :u)
channels(V::OneDimension) = (:τ,)
channels(Vs::FourFields) = (:s, :t, :u)
channels(V::OneField) = (:τ,)

function CorrelationNonChiral(V::ExtFields{T}, Nmax::Int) where {T}
    @assert all((v.c === V[1].c for v in V)) """
    External fields in the argument of the Correlation constructor do not all have the same
    CentralCharge
    """

    dims = Tuple(
        Tuple(v.dims[lr] for v in V)
        for lr in (:left, :right)
    )
    corr_chiral = Tuple(
        CorrelationChiral(dims[lr], Nmax)
        for lr in (:left, :right)
    )

    Rmn = Tuple(corr_chiral[lr]._Rmn for lr in (:left, :right))
    Rmnreg = Tuple(corr_chiral[lr]._Rmn_reg for lr in (:left, :right))
    CNmn = Tuple(corr_chiral[lr]._CNmn for lr in (:left, :right))

    CorrelationNonChiral{T}(V, Nmax, Rmn, Rmnreg, CNmn)
end

function Base.getproperty(co::CorrelationNonChiral, s::Symbol)
    s === :c && return getfield(co, :fields)[1].c
    return getfield(co, s)
end

function Base.getindex(c::CorrelationNonChiral, s::Symbol)
    s in (:left, :right) && return CorrelationChiral(
        Tuple(v.dims[s] for v in c.fields),
        c.Nmax, c._Rmn[s], c._Rmn_reg[s], c._CNmn[s]
    )
    error("cannot access Correlation")
end

function CorrelationNonChiral(corr_chiral::LeftRight{CorrelationChiral{T}}) where {T}
    dims_left = corr_chiral[:left].dims
    dims_right = corr_chiral[:left].dims
    fields = Tuple(
        Field(d1, d2)
        for (d1, d2) in zip(dims_left, dims_right)
    )
    Nmax = corr_chiral[:left].Nmax
    Rmn = Tuple(corr_chiral[lr]._Rmn for lr in (:left, :right))
    Rmnreg = Tuple(corr_chiral[lr]._Rmn_reg for lr in (:left, :right))
    CNmn = Tuple(corr_chiral[lr]._CNmn for lr in (:left, :right)) 

    CorrelationNonChiral(fields, Nmax, Rmn, Rmnreg, CNmn)
end

function Correlation()
    CorrelationNonChiral(Field(), 0)
end

Correlation(ds::ExtDimensions, Nmax::Int) = CorrelationChiral(ds, Nmax)

Correlation(Vs::ExtFields, Nmax::Int) = CorrelationNonChiral(Vs, Nmax)

Correlation(Vs::ExtFields, lr::Symbol, Nmax::Int) = CorrelationChiral(
    Tuple(v.dims[lr] for v in Vs), Nmax
)

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

Correlation(args...; Δmax=10.) = Correlation(args..., N_max(CD(), Δmax))

