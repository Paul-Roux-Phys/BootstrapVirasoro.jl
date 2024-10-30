const RmnTable{T}  = Dict{Tuple{Int, Int},      T}
const CNmnTable{T} = Dict{Tuple{Int, Int, Int}, T}

const Channels{T} = Dict{Symbol, T} # type for holding data in all channels

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
    _Rmn::Channels{LeftRight{RmnTable{T}}}
    _Rmn_reg::Channels{LeftRight{RmnTable{T}}}
    _CNmn::Channels{LeftRight{CNmnTable{T}}}

end

function Base.show(io::IO, corr::Correlation)
    print("Correlation function with external fields\n$(corr.fields)")
end

function permute_fields(Vs::FourFields, chan::Symbol)
    chan === :s && return Vs
    chan === :t && return (Vs[1], Vs[4], Vs[3], Vs[2])
    chan === :u && return (Vs[1], Vs[3], Vs[2], Vs[4])
end

permute_fields(V::OneField, chan::Symbol) = V

function channels(V::FourFields)
    return (:s, :t, :u)
end

function channels(V::OneField)
    return (:τ,)
end

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
    Rmn = Channels{LeftRight{RmnTable{T}}}()
    Rmnreg = Channels{LeftRight{RmnTable{T}}}()
    CNmn = Channels{LeftRight{CNmnTable{T}}}()

    for x in channels(V)

        Vx = permute_fields(V, x)

        Rmn[x] = Tuple(
            Dict(
                (m, n) => computeRmn(m, n, Vx, lr)
                for m in 1:Nmax, n in 1:Nmax
                if Rmn_zero_order(m, n, Vx) == 0 && m * n <= Nmax
            )
            for lr in (:left, :right)
        )

        # Rmnreg[x] = LeftRight(RmnTable{T}() for lr in (:left, :right))

        Rmnreg[x] = Tuple(
            Dict(
                (m, n) => computeRmn(m, n, Vx, lr)
                for m in 1:Nmax, n in 1:Nmax
                if Rmn_zero_order(m, n, Vx) > 0 && m*n <= Nmax
            )
            for lr in (:left, :right)
        )

        CNmn[x] = Tuple(
            Dict(
                (N, m, n) => computeCNmn(N, m, n, V[1].c, Rmn[x][lr])
                for m in 1:Nmax, n in 1:Nmax, N in 1:Nmax
                if m * n <= N && (m, n) in keys(Rmn[x][lr])
            )
            for lr in (:left, :right)
        )
        
    end

    Correlation{T}(V, Nmax, Rmn, Rmnreg, CNmn)
end

function Base.getproperty(co::Correlation, s::Symbol)
    s === :c && return co.fields[1].c
    return getfield(co, s)
end

Correlation(V1, V2, V3, V4, Nmax::Int) = Correlation((V1, V2, V3, V4), Nmax)
Correlation(V, Nmax) = Correlation((V,), Nmax)