"""
    Field{T}
Object representing a conformal field.
Contains the conformal dimensions, and flags saying whether the field has rational Kac
indices and/or is diagonal. The user can access different parametrisations of the left and 
right conformal dimensions, and Kac indices.

# Examples

```jldoctest
julia> c = CentralCharge(B = 0.5)
c = 27.999999999999996 + 0.0im, β = 0.0 - 0.7071067811865476im

julia> V = Field(c, r=0, s=2//3)
Diagonal Field{ComplexF64} with ConformalDimension{ComplexF64} with Kac indices r = 0//1, s = 2//3

julia> V.r
0//1

julia> V.s
2//3

julia> V.δ
(-0.22222222222222215 - 0.0im, -0.22222222222222215 + 0.0im)

julia> V.δ[:left]
-0.22222222222222215 - 0.0im

julia> V.p[:right]
0.4714045207910316 + 0.0im

julia> V2 = Field(c, δ = 0.5, diagonal=true)
Diagonal Field{ComplexF64} with ConformalDimension{ComplexF64} with
Δ = 1.625 + 0.0im, P = 0.7071067811865476
```
"""
struct Field{T}

    dims::LeftRight{ConformalDimension{T}}

end

"""
```julia
    Field(charge, parameter = value; kwargs)
    Field(charge, r = r, s = s)
    Field((dim_left, dim_right))
    Field(dim_left, dim_right)
    Field(dim) # diagonal field
```

Constructor function for the Field type.

Given a charge `charge`, one of the four parameters `Δ`, `δ`, `P`, `p` and two values,
create an object `Field{T}` (where T is the type of the values in `charge`) that represents
a field of left and right dimensions given by leftvalue and rightvalue in the chosen
parametrisation.

# keyword arguments:

- `r::Rational`,`s::Rational`: used in conjunction to `Kac=true`, must be given rational
values. By convention ``V_(r,s)`` has left and right momenta ``(P_(r,s), P_(r,-s))``.
- `diagonal::Bool`: set to `true` to get a diagonal field;

# Examples
```jldoctest
julia> setprecision(BigFloat, 20, base=10);

julia> c = CentralCharge(β = big"0.5");

julia> V = Field(c, r=0, s=1)
Diagonal Field{Complex{BigFloat}} with ConformalDimension{Complex{BigFloat}} with Kac indices r = 0, s = 1

julia> V.Δ
(0.4375 + 0.0im, 0.4375 + 0.0im)

julia> V.P[:left]
1.0

julia> V.p[:right]
-0.0 + 1.0im

julia> V2 = Field(c, :P, 0.42, diagonal=true)
Diagonal Field{Complex{BigFloat}} with ConformalDimension{Complex{BigFloat}} with
Δ = -0.3861000000000000130545 + 0.0im, P = 0.4199999999999999844569
```
"""
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
    if diagonal
        dim_right = dim_left
    else
        @assert (r !== missing && s !== missing) "
            A non-diagonal field must be given from Kac indices.
            If you mean to define a diagonal field, use `diagonal=true`.
        "
        r*s % 1 != 0 && @info "You defined a field with non-integer r*s, is that intentional?"
        dim_right = ConformalDimension(c, sym, dim_left, r=r, s=-s)
    end

    Field{T}(LeftRight((dim_left, dim_right)))
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
function Field(ds::LeftRight{ConformalDimension})
    diagonal = false
    degenerate = false
    if ds[:left] == ds[:right]
        diagonal = true
    end
    if ds[:left].isKac && ds[:right].isKac 
        degenerate=true
    end
    return Field(ds, degenerate)
end
Field(d_left::ConformalDimension, d_right::ConformalDimension) = Field((d_left, d_right))
Field(d::ConformalDimension) = Field(d, d)

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
    s === :dim && getproperty(V, :r) == 0 && return ds[:left]

    return getfield(V, s)
end

# Overload the == operator
function Base.:(==)(V1::Field, V2::Field)
    return V1.Δ == V2.Δ
end

"""Compute the spin Δleft - Δright of a field."""
function spin(V::Field)::Rational
    if V.r == 0
        return 0
    elseif V.isKac
        return V.r*V.s
    else # this should never happen
        @warn "You are computing the spin of a field not defined from Kac indices"
        return V.Δ[1] - V.Δ[2]
    end
end

function swap_lr(V::Field{T}) where {T}
    return Field{T}((V.dims[:right], V.dims[:left]))
end

function Base.show(io::IO, ::MIME"text/plain", V::Field)
    if V.r == 0
        print(io, "Diagonal $(typeof(V)) with ")
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
    if V.r == 0
        print(io, "V_{P=$(V.P[:left])}")
    else
        print(io, "V_{$(V.indices)}")
    end
end

function shift(V::Field, i, index=:s)
    if V.r == 0
        Field(shift(V.dim, i, index))
    else
        Field(shift(V.dims[:left], i, index), shift(V.dims[:right], -i, index))
    end
end

scaling_dim(V::Field) = V.Δ[:left] + V.Δ[:right]
