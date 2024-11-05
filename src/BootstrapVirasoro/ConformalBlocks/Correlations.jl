const RmnTable{T}  = Dict{Tuple{Int, Int},      T}
const CNmnTable{T} = Dict{Tuple{Int, Int, Int}, T}

const Channels{T} = Dict{Symbol, T} # type for holding data in all channels

struct CorrelationChiral{T}

    dims::ExtDimensions{T}
    Nmax::Int
    _Rmn::Channels{RmnTable{T}}
    _Rmn_reg::Channels{RmnTable{T}}
    _CNmn::Channels{CNmnTable{T}}

end

function Base.show(io::IO, corr::CorrelationChiral)
    print("Chiral correlation function with external dimensions\n$(corr.dims)")
end

function CorrelationChiral(d::ExtDimensions{T}, Nmax::Int) where {T}
    @assert all((dim.c === d[1].c for dim in d)) """
    External fields in the argument of the Correlation constructor do not all have the same
    CentralCharge
    """
    Rmn = Channels{RmnTable{T}}()
    Rmnreg = Channels{RmnTable{T}}()
    CNmn = Channels{CNmnTable{T}}()

    for x in channels(d)

        dx = permute_dimensions(d, x)

        Rmn[x] = Dict(
            (m, n) => computeRmn(m, n, dx)
            for m in 1:Nmax, n in 1:Nmax
            if Rmn_zero_order(m, n, dx) == 0 && m * n <= Nmax
        )

        Rmnreg[x] = Dict(
            (m, n) => computeRmn(m, n, dx)
            for m in 1:Nmax, n in 1:Nmax
            if Rmn_zero_order(m, n, dx) > 0 && m * n <= Nmax
        )

        CNmn[x] = Dict(
            (N, m, n) => computeCNmn(N, m, n, d[1].c, Rmn[x])
            for m in 1:Nmax, n in 1:Nmax, N in 1:Nmax
            if m * n <= N && (m, n) in keys(Rmn[x])
        )

    end

    CorrelationChiral{T}(d, Nmax, Rmn, Rmnreg, CNmn)
end

CorrelationChiral(d1, d2, d3, d4, Nmax) = CorrelationChiral((d1, d2, d3, d4), Nmax)

"""
Correlation{T}

Struct for holding data common to all conformal blocks appearing in a given correlation
function:
* External fields
* Residues `R_mn` in all channels
* Regularised residues `R_mn` when some R_mn are zero
* Coefficients `C^N_{mn}` in all channels
"""
struct Correlation{T}

    fields::ExtFields{T}
    Nmax::Int
    _Rmn::LeftRight{Channels{RmnTable{T}}}
    _Rmn_reg::LeftRight{Channels{RmnTable{T}}}
    _CNmn::LeftRight{Channels{CNmnTable{T}}}

end

function Base.show(io::IO, corr::Correlation)
    print("Correlation function with external fields\n$(corr.fields)")
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

"""
    Correlation(Vs::FourFields{T}, Nmax::Int)

Constructor function for the Correlation type.

# Examples

```julia-repl
julia> c = CentralCharge(:c, 0.5);
julia> V1 = Field(c, Kac=true, r=1, s=0);
julia> V2 = Field(c, Kac=true, r=2, s=0);
julia> corr = Correlation((V1, V2, V2, V1), 6)
```

Note that the parentheses are optional around the four fields:
```julia-repl
julia> corr2 = Correlation(V1, V1, V2, V1, 6);
julia> corr2._Rmn[:u][:left]
```
"""
function Correlation(V::ExtFields{T}, Nmax::Int) where {T}
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

    Correlation{T}(V, Nmax, Rmn, Rmnreg, CNmn)
end

function Base.getproperty(co::Correlation, s::Symbol)
    s === :c && return getfield(co, :fields)[1].c
    return getfield(co, s)
end

function Base.getindex(c::Correlation, s::Symbol)
    s in (:left, :right) && return CorrelationChiral(
        Tuple(v.dims[s] for v in c.fields),
        c.Nmax, c._Rmn[s], c._Rmn_reg[s], c._CNmn[s]
    )
end

function Correlation(corr_chiral::LeftRight{CorrelationChiral{T}}) where {T}
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

    Correlation(fields, Nmax, Rmn, Rmnreg, CNmn)
end
Correlation(co_left::CorrelationChiral, co_right::CorrelationChiral) = Correlation((co_left, co_right))
Correlation(V1, V2, V3, V4, Nmax::Int) = Correlation((V1, V2, V3, V4), Nmax)
Correlation(V, Nmax) = Correlation((V,), Nmax)