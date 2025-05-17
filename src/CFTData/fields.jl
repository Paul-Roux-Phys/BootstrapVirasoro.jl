"""
    Field(c, parameter=value; r, s, diagonal=false)

Type for representing a non-chiral field.

# keyword arguments:

- `r::Rational`,`s::Rational`. By convention ``V_{(r,s)}`` has left and right momenta ``(P_{(r,s)}, P_{(r,-s)})``.
- `diagonal::Bool`: set to `true` to get a diagonal field;

# Examples

```jldoctest
julia> c = CentralCharge(β = big"0.5");

julia> V = Field(c, r=0, s=1)
Diagonal Field{Complex{BigFloat}}, dim = Δ_{0, 1}

julia> V.Δ
(0.4375 + 0.0im, 0.4375 + 0.0im)

julia> V.P[:left] ≈ 1
true

julia> V.p[:right] ≈ -im
true

julia> V2 = Field(c, :P, 0.42, diagonal=true); isdiagonal(V2)
true
```
"""
struct Field{T}

    dims::LeftRight{ConformalDimension{T}}
    diagonal::Bool

end

function Field(
    c::CentralCharge{T},
    sym::Symbol,
    dim;
    r=missing, s=missing,
    degenerate=false, diagonal=false
) where {T}
    if degenerate # degenerate fields are diagonal and must be given from Kac indices
        @assert (r !== missing && s !== missing) "
            A degenerate field must be given from Kac indices.
        "
        diagonal = true
    end
    dim_left = ConformalDimension(c, sym, dim, r=r, s=s)
    if diagonal || r == 0
        diagonal = true
        dim_right = dim_left
    else
        @assert (r !== missing && s !== missing) "
            A non-diagonal field must be given from Kac indices.
            If you mean to define a diagonal field, use `diagonal=true`.
        "
        r*s % 1 != 0 && @info "You defined a field with non-integer r*s, is that intentional?"
        dim_right = ConformalDimension(c, sym, dim_left, r=r, s=-s)
    end

    Field{T}(LeftRight((dim_left, dim_right)), diagonal)
end

function Field(
    c::CentralCharge;
    r=missing, s=missing,
    diagonal=false, degenerate=false,
    Δ=missing, δ=missing, P=missing, p=missing
)
    r !== missing && s !== missing &&
        return Field(c, :Δ, 0, r=r, s=s, degenerate=degenerate, diagonal=diagonal)
    Δ !== missing && return Field(c, :Δ, Δ, diagonal=true)
    δ !== missing && return Field(c, :δ, δ, diagonal=true)
    P !== missing && return Field(c, :P, P, diagonal=true)
    p !== missing && return Field(c, :p, p, diagonal=true)
    return Field(c, :Δ, 0, r=1, s=1, diagonal=true)
end

Field() = Field(CentralCharge())
function Field(ds::LeftRight{ConformalDimension}; diagonal=false)
    return Field(ds, diagonal)
end
Field(d_left::ConformalDimension, d_right::ConformalDimension; diagonal=false) = Field((d_left, d_right); diagonal=diagonal)
Field(d::ConformalDimension; diagonal=true) = Field(d, d, diagonal=diagonal)

function Base.getproperty(V::Field, s::Symbol)
    ds = getfield(V, :dims)
    s === :c && return V.dims[:left].c
    s === :P && return LeftRight((ds[:left].P, ds[:right].P))
    s === :Δ && return LeftRight((ds[:left].Δ, ds[:right].Δ))
    s === :p && return LeftRight((ds[:left].p, ds[:right].p))
    s === :δ && return LeftRight((ds[:left].δ, ds[:right].δ))
    (s === :r || s === :s) && return getproperty(ds[:left], s) # by convention V_(r,s) denotes the field
                                             # with left right dimension P_(r, s), P_(r, -s)
    s === :isKac && return (V.dims[:left].isKac && V.dims[:right].isKac)
    s === :indices && return ds[:left].indices
    s === :dim && isdiagonal(V) && return ds[:left]

    return getfield(V, s)
end

# Overload the == operator
function Base.:(==)(V1::Field, V2::Field)
    return V1.Δ == V2.Δ
end

"""
        isdiagonal(V)::Bool
Whether V is a diagonal field.
"""
isdiagonal(V::Field) = V.diagonal

"""
        isdegenerate(V)::Bool
Whether V is a degenerate field.
"""
function isdegenerate(V::Field)
    return V.r % 1 == 0 && real(V.s) % 1 == 0 &&
        imag(V.s) == 0 &&
        V.r > 0 && V.s > 0 && isdiagonal(V) 
end

"""Spin(V::Field) = Δleft - Δright."""
function spin(V::Field)::Rational
    if isdiagonal(V)
        return 0
    elseif V.isKac
        return V.r*V.s
    else # this should never happen
        @warn "You are computing the spin of a field not defined from Kac indices"
        return V.Δ[1] - V.Δ[2]
    end
end

"""
        swap_lr(V)
Return a field with left and right dimensions swapped.
"""
function swap_lr(V::Field{T}) where {T}
    return Field{T}((V.dims[:right], V.dims[:left]))
end

"""
        shift(V, i)

Shift the field:
- s -> s+i if !isdiagonal(V)
- P -> P+i/(2β) if isdiagonal(V)
"""
function shift(V::Field, i, index=:s)
    if isdiagonal(V)
        Field(shift(V.dim, i, index), diagonal=true)
    else
        Field(shift(V.dims[:left], i, index), shift(V.dims[:right], -i, index))
    end
end

"""
        total_dimension(V)
Return ``Δ + \barΔ``.
"""
total_dimension(V::Field) = V.Δ[:left] + V.Δ[:right]

function Base.hash(V::Field, h::UInt)
    return hash(V.dims, h)
end

function Base.show(io::IO, ::MIME"text/plain", V::Field)
    if isdiagonal(V)
        print(io, "Diagonal $(typeof(V)), dim = ")
        show(io, V.dims[:left])
    else
        print(io, "Non-diagonal $(typeof(V))\n")
        print(io, "left: ")
        show(io, V.dims[:left])
        print(io, "\nright: ")
        show(io, V.dims[:right])
    end
end

function Base.show(io::IO, V::Field)
    if isdiagonal(V)
        if V.isKac
            if isdegenerate(V)
                print(io, "<$(V.r), $(V.s)>")
            else
                print(io, "($(V.r), $(V.s))")
            end
        else
            print(io, "V_{P=$(V.P[:left])}")
        end
    else
        print(io, "V_{$(V.indices)}")
    end
end
