"""
    Field(c, parameter=value; r, s, diagonal=false)

Type for representing a non-chiral field.

# keyword arguments:

- `r::Rational`,`s::Rational`. By convention ``V_{(r,s)}`` has
left and right momenta ``(P_{(r,s)}, P_{(r,-s)})``.
- `diagonal::Bool`: set to `true` to get a diagonal field;

# Examples

```jldoctest
julia> c = CentralCharge(β = big"0.5");

julia> V = Field(c, r=0, s=1)
Diagonal Field{Complex{BigFloat}}, dim = Δ_{0, 1}

julia> (V.dims[:left].Δ, V.dims[:right].Δ)
(0.4375 + 0.0im, 0.4375 + 0.0im)

julia> V.dims[:left].P ≈ 1
true

julia> V.dims[:right].p ≈ -im
true

julia> V2 = Field(c, :P, 0.42, diagonal=true); isdiagonal(V2)
true
```
"""
struct Field{T}
    c::CentralCharge{T}
    dims::LeftRight{ConformalDimension{T}}
    r::Union{Rational,Int}
    s::Union{T,Rational,Int}
    diagonal::Bool
    degenerate::Bool
    isKac::Bool
end

function Field(
    c::CC{T},
    sym::Symbol,
    dim;
    r = missing,
    s = missing,
    degenerate = false,
    diagonal = false,
) where {T}
    if degenerate # degenerate fields are diagonal and must be given from Kac indices
        @assert (r !== missing && s !== missing) "
            A degenerate field must be given from Kac indices.
        "
        diagonal = true
    end
    if ismissing(r) ⊻ ismissing(s) # xor
        error("You cannot give only r or only s, you must give both")
    end
    if ismissing(r) && ismissing(s) || r == 0
        diagonal = true
    end
    dim_left = ConformalDimension(c, sym, dim, r = r, s = s)
    if diagonal
        diagonal = true
        dim_right = dim_left
    else
        @assert (r !== missing && s !== missing) "
            A non-diagonal field must be given from Kac indices.
            If you mean to define a diagonal field, use `diagonal=true`.
        "
        r * s % 1 != 0 &&
            @info "You defined a field with non-integer r*s, is that intentional?"
        dim_right = ConformalDimension(c, sym, dim_left, r = r, s = -s)
    end
    if ismissing(r) || ismissing(s)
        r, s = 0, 2c.β * dim_left.P
    end
    degenerate =
        r % 1 == 0 && real(s) % 1 == 0 && imag(s) == 0 && r > 0 && s > 0 && diagonal
    r, s = indices(dim_left)
    isKac = dim_left.isKac && dim_right.isKac

    Field{T}(c, LeftRight((dim_left, dim_right)), r, s, diagonal, degenerate, isKac)
end

function Field(
    c::CC;
    r = missing,
    s = missing,
    diagonal = false,
    degenerate = false,
    Δ = missing,
    δ = missing,
    P = missing,
    p = missing,
)
    r !== missing &&
        s !== missing &&
        return Field(
            c,
            :Δ,
            0,
            r = r,
            s = s,
            degenerate = degenerate,
            diagonal = diagonal,
        )
    Δ !== missing && return Field(c, :Δ, Δ, diagonal = true)
    δ !== missing && return Field(c, :δ, δ, diagonal = true)
    P !== missing && return Field(c, :P, P, diagonal = true)
    p !== missing && return Field(c, :p, p, diagonal = true)
    return Field(c, :Δ, 0, r = 1, s = 1, diagonal = true)
end

Field() = Field(CentralCharge())
function Field(ds::LeftRight{CD{T}}) where {T}
    c = ds[1].c
    degenerate = isdegenerate(ds[1]) && isdegenerate(ds[2])
    isKac = ds[1].isKac && ds[2].isKac
    ds[1] == ds[2] ? diagonal=true : diagonal = false;
    r, s = ds[1].r, ds[1].s
    Field{T}(c, ds, r, s, diagonal, degenerate, isKac)
end
Field(d1::CD, d2::CD) = Field((d1, d2)::LeftRight)
Field(d::CD) = Field((d, d)::LeftRight)

"""
        isdiagonal(V)::Bool
Whether V is a diagonal field.
"""
isdiagonal(V::Field) = V.diagonal

"""
        isdegenerate(V)::Bool
Whether V is a degenerate field.
"""
isdegenerate(V::Field) = V.degenerate

indices(V::Field) = V.r, V.s

"""Spin(V::Field) = Δleft - Δright."""
function spin(V::Field)::Int
    if isdiagonal(V)
        return 0
    elseif V.isKac
        return Int(V.r * V.s)
    else # this should never happen
        @warn "You are computing the spin of a field not defined from Kac indices"
        return V.dims[1].Δ - V.dims[2].Δ
    end
end

"""
        swap_lr(V)
Return a field with left and right dimensions swapped.
"""
function swap_lr(V::Field{T}) where {T}
    return Field((V.dims[:right], V.dims[:left]))
end

"""
        shift(V, i)

Shift the field:
- s -> s+i if !isdiagonal(V)
- P -> P+i/(2β) if isdiagonal(V)
"""
function shift(V::Field, i, index = :s)
    if isdiagonal(V)
        Field(shift(V.dims[:left], i, index))
    else
        Field(shift(V.dims[:left], i, index), shift(V.dims[:right], -i, index))
    end
end

reflect(V) = Field(V.c, r = V.r, s = -V.s)

"""
        total_dimension(V)
Return ``Δ + \barΔ``.
"""
total_dimension(V::Field) = V.dims[:left].Δ + V.dims[:right].Δ

# Implement the hashing interface (for Dict, Set)
function Base.isequal(a::Field, b::Field)
    d = (a.diagonal == b.diagonal)
    l = isequal(a.dims[:left], b.dims[:left])
    r = isequal(a.dims[:right], b.dims[:right])
    return d && l && r
end

Base.:(==)(V1::Field, V2::Field) = isequal(V1, V2)

function Base.hash(V::Field, h::UInt)
    # hash using the string representation
    str = sprint(show, V)
    return hash(str, h)
end

function Base.show(io::IO, ::MIME"text/plain", V::Field)
    if isdegenerate(V)
        print(io, "Degenerate $(typeof(V)), dim = ")
        show(io, V.dims[:left])
    elseif isdiagonal(V)
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
                print(io, "V_{<$(V.r), $(V.s)>}")
            else
                print(io, "V_{($(V.r), $(V.s))}")
            end
        else
            print(io, "V_{P=$(V.dims[:left].P)}")
        end
    else
        print(io, "V_{$(indices(V))}")
    end
end
