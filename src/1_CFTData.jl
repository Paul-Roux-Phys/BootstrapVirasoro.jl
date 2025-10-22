"""
        CentralCharge(; β, c, B, b)

Create a CentralCharge object from one of the parameters 
``β, c = 13 - 6 β^2 - 6 β^{-2}, b = i β, B = b^2``.

Aliased to CC.

# Examples

```jldoctest; output=false
c1 = CentralCharge(β=0.5+0.3im)
c2 = CC(c=big"0.4")
c2.c ≈ big"0.4" # access any of the parameters

# output

true
```
"""
struct CentralCharge{T}
    β::T
    B::T
    b::T
    c::T
    n::T
end
# convenience alias
const CC = CentralCharge


"""
        ConformalDimension(c::CC; Δ, δ, P, p, r, s)

Create a ConformalDimension object from one of the parameters 
``P, p = iP, δ = P^2,  Δ = δ + \\frac{c-1}{24}`` or from Kac indices
`r`, `s` (``P_{(r, s)} = \\frac{1}{2}(rβ - sβ^{-1})``). The index `r` must be integer or rational, `s` can be arbitrary.

Aliased to CD.

# Examples

```jldoctest; output=false
c = CC(β = big"0.3"+big"0.1"*im)
d1 = CD(c, P=big"0.5"*im)
d2 = CD(c, r=2, s=1//2)
d3 = CD(c, r=0, s=big"0.1")
d2.r == 2 # true

# output

true
```
"""
struct ConformalDimension{T}
    c::CentralCharge{T}
    P::T
    p::T
    δ::T
    Δ::T
    r::Union{Rational,Int}   # r must be rational
    s::Union{T,Rational,Int} # s can be an arbitrary number to represent a
    isKac::Bool              # diagonal field. Still needs to be exact in
    degenerate::Bool         # case of a non-diagonal field
end
# convenience alias
const CD = ConformalDimension

"""
        Field(c::CC; Δ, δ, P, p, r, s, diagonal)

Create a Field object from one of the parameters 
`P, p, δ, Δ` (see [ConformalDimension](@ref)) or from Kac indices
`r`, `s`. The index `r` must be integer or rational, `s` can be arbitrary.
Left and right dimensions are accessed with [:left], [:right].

# Examples

```jldoctest; output=false
c = CC(β = big"0.3"+big"0.1"*im)
V1 = Field(c, r=2, s=1//2) # non-diagonal, (Δ(2, 1//2), Δ(2, -1//2))
V2 = Field(c, P=0.5)       # diagonal
V3 = Field(c, r=1, s=2, diagonal=true) # degenerate
V3[:left].Δ == V3[:right].Δ # true

# output

true
```
"""
struct Field{T}
    c::CC{T}
    dims::LeftRight{ConformalDimension{T}}
    r::Union{Rational,Int}
    s::Union{T,Rational,Int}
    diagonal::Bool
    degenerate::Bool
    isKac::Bool
end

#======================================================================================
Central charges
======================================================================================#
function Bfrom(s::Symbol, x)
    a = (x-1)*(x-25)
    s === :β && return -x^2
    s === :c && return (
        if a isa Real && a > 0
            (x - 13 + sqrt((x - 1) * (x - 25))) / 12
        else # a is complex
            (x - 13 + sqrt(complex((x - 1) * (x - 25)))) / 12
        end
    )
    s === :b && return x^2
    s === :B && return x
    error("unsupported parameter: $s")
end

function Bto(s::Symbol, x)
    rx = sqrt(complex(x))
    s === :β && return im*rx
    s === :c && return 13+6*x+6/x
    s === :b && return rx
    s === :B && return x
    s === :n && return -2cos(π*x)
end

function _CentralCharge(s::Symbol, x)
    B = complex(Bfrom(s, x))
    β = Bto(:β, B)
    b = Bto(:b, B)
    c = Bto(:c, B)
    n = Bto(:n, B)
    CentralCharge(β, B, b, c, n)
end

function CentralCharge(;
    β = missing,
    c = missing,
    B = missing,
    b = missing,
    n = missing,
)
    β !== missing && return _CentralCharge(:β, β)
    c !== missing && return _CentralCharge(:c, c)
    B !== missing && return _CentralCharge(:B, B)
    b !== missing && return _CentralCharge(:b, b)
    n !== missing && return _CentralCharge(:n, n)
    return _CentralCharge(:c, 1)
end

function Base.show(io::IO, ::MIME"text/plain", c::CC)
    print(io, "c = $(c.c), β = $(c.β)")
end

function Base.show(io::IO, c::CC)
    print(io, "c = $(c.c)")
end

function Base.isequal(c1::CC, c2::CC)
    return c1.β == c2.β
end

function Base.hash(c::CC, h::UInt)
    return hash(c.β, h)
end

#======================================================================================
Conformal Dimensions
======================================================================================#
function Pfrom(s::Symbol, x, c::CC)
    s === :Δ && return sqrt(complex(x - (c.c-1)/24))
    s === :δ && return sqrt(complex(x))
    s === :P && return x
    s === :p && return im*x
end

function Pto(s::Symbol, x, c::CC)
    s === :Δ && return x^2 + (c.c-1)/24
    s === :δ && return x^2
    s === :P && return x
    s === :p && return -im*x
    s === :w && return -2*cos(oftype(c.β, π)*c.β*x)
end

Prs(r, s, β::Number) = (β*r - s/β) / 2
"""
``P_{(r, s)} = 1/2 * (β r - s / β)``.
"""
Prs(r, s, c::CC) = Prs(r, s, c.β)
δrs(r, s, B::Number) = -1/4 * (B*r^2 + 2*r*s + s^2/B)
δrs(r, s, c::CC) = -1/4 * (c.B*r^2 + 2*r*s + s^2/c.B)

function _ConformalDimension(c::CC{T}, sym::Symbol, P, r, s) where {T<:Number}
    β = c.β
    degenerate = false
    isKac = false
    if (r !== missing && s !== missing)
        r isa Real && r % 1 == 0 ? r = Int(r) : nothing;
        s isa Real && s % 1 == 0 ? s = Int(s) : nothing;
        P = Prs(r, s, β)
        if (r isa Rational || r isa Integer) && (s isa Rational || s isa Integer)
            isKac = true
            if (r isa Integer && r > 0 && s isa Integer && s > 0)
                degenerate = true
            end
        end
        if r == 0 && !(s isa Int || s isa Rational)
            s = convert(T, s)
        end
    else
        @assert (r === missing && s === missing) "
            You cannot give only r or only s, you must give neither or both
        "
        P = Pto(:P, Pfrom(sym, P, c), c)
        r = 0
        s = - 2β * P
    end
    p = Pto(:p, P, c)
    δ = Pto(:δ, P, c)
    Δ = Pto(:Δ, P, c)
    ConformalDimension{T}(c, P, p, δ, Δ, r, s, isKac, degenerate)
end

function ConformalDimension(
    c::CC;
    r = missing,
    s = missing,
    Δ = missing,
    δ = missing,
    P = missing,
    p = missing,
)
    (r !== missing && s !== missing) && return _ConformalDimension(c, :Δ, 0, r, s)
    Δ !== missing && return _ConformalDimension(c, :Δ, Δ, r, s)
    δ !== missing && return _ConformalDimension(c, :δ, δ, r, s)
    P !== missing && return _ConformalDimension(c, :P, P, r, s)
    p !== missing && return _ConformalDimension(c, :p, p, r, s)
    return _ConformalDimension(c, :Δ, 0, r, s)
end

ConformalDimension() = ConformalDimension(CC())

indices(d::CD) = d.r, d.s

function Base.:+(d1::CD, d2::CD)
    ConformalDimension(d1.c, Δ = d1.Δ + d2.Δ)
end

""" s → s+shift, P->P-shift/2/β """
shift(d::CD, shift) = ConformalDimension(d.c, r = d.r, s = d.s + shift)

function Base.isequal(a::CD, b::CD)
    c = isequal(a.c, b.c)
    p = isequal(a.P, b.P)
    k = (a.isKac == b.isKac) && a.r == b.r && a.s == b.s
    return c && p && k
end

Base.:(==)(V1::CD, V2::CD) = isequal(V1, V2)

function Base.hash(d::CD, h::UInt)
    str = sprint(show, d)
    return hash(str, h)
end

function Base.show(io::IO, ::MIME"text/plain", d::CD{T}) where {T}
    if d.isKac
        print(io, "ConformalDimension{$T} with Kac indices r = $(d.r), s = $(d.s)")
    else
        print(io, "ConformalDimension{$T} with\nΔ = $(d.Δ), P = $(d.P)")
    end
end

function Base.show(io::IO, d::CD{T}) where {T}
    if d.isKac
        print(io, "Δ_{$(d.r), $(d.s)}")
    else
        print(io, "Δ_{P=$(d.P)}")
    end
end

#======================================================================================
Fields
======================================================================================#
function _Field(c::CC{T}, sym::Symbol, dim, r, s, diagonal) where {T}
    @assert !(ismissing(r) ⊻ ismissing(s)) "
    You cannot give only r or only s, you must give neither or both
"

    if ismissing(r) || r == 0
        diagonal = true
    end
    dim_left = _ConformalDimension(c, sym, dim, r, s)

    if diagonal
        dim_right = dim_left
    else
        @assert (r !== missing && s !== missing) "
            A non-diagonal field must be given from Kac indices.
            If you mean to define a diagonal field, use `diagonal=true`.
        "
        dim_right = ConformalDimension(c, r = r, s = -s)
    end
    if ismissing(r)
        r, s = 0, -2c.β * dim_left.P
    end
    r, s = indices(dim_left)
    isKac = dim_left.isKac && dim_right.isKac
    degenerate = dim_left.degenerate && dim_right.degenerate

    Field{T}(
        c,
        LeftRight{CD{T}}(dim_left, dim_right),
        r,
        s,
        diagonal,
        degenerate,
        isKac,
    )
end

function Field(
    c::CC;
    r = missing,
    s = missing,
    diagonal = false,
    Δ = missing,
    δ = missing,
    P = missing,
    p = missing,
)
    r !== missing && s !== missing && return _Field(c, :Δ, 0, r, s, diagonal)
    Δ !== missing && return _Field(c, :Δ, Δ, missing, missing, true)
    δ !== missing && return _Field(c, :δ, δ, missing, missing, true)
    P !== missing && return _Field(c, :P, P, missing, missing, true)
    p !== missing && return _Field(c, :p, p, missing, missing, true)
    return _Field(c, :Δ, 0, 1, 1, true)
end

function Field(ds::LeftRight{CD{T}}) where {T}
    c = ds.left.c
    degenerate = ds.left.degenerate && ds.right.degenerate
    isKac = ds.left.isKac && ds.right.isKac
    ds.left == ds.right ? diagonal=true : diagonal = false;
    r, s = indices(ds.left)
    Field{T}(c, ds, r, s, diagonal, degenerate, isKac)
end
Field(ds::Tuple{CD{T},CD{T}}) where {T} = Field(LeftRight{CD{T}}(ds[1], ds[2]))
Field(d1::CD{T}, d2::CD) where {T} = Field(LeftRight{CD{T}}(d1, d2))
Field(d::CD{T}) where {T} = Field(LeftRight{CD{T}}(d, d))
Field() = Field(CC())

Base.getindex(V::Field, s::Symbol) = getfield(V.dims, s)

indices(V::Field) = V.r, V.s
function get_indices(V)
    # not the same as indices because it returns (0, -2βP) also for degenerate fields
    # the latter returns integer (r, s) for the degenerate fields.
    β = V.c.β
    if V.diagonal
        return 0, -2*β*V[:left].P
    else
        return indices(V)
    end
end

"""Spin(V::Field) = Δleft - Δright."""
function spin(V::Field{T})::Union{Int,Rational,T} where {T}
    V.diagonal && return 0
    spin = V.isKac ? V.r * V.s : V[:left].Δ - V.dims.right.Δ
    return spin % 1 == 0 ? Int(spin) : spin
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
- s -> s+i if !V.diagonal
- P -> P+i/(2β) if V.diagonal
"""
function shift(V::Field, i)
    if V.diagonal
        Field(shift(V.dims[:left], i))
    else
        Field(shift(V.dims[:left], i), shift(V.dims[:right], -i))
    end
end

function reflect(V)
    V.diagonal && return V
    return Field(V.c, r = V.r, s = -V.s)
end

""" ``Δ + \barΔ``. """
total_dimension(V::Field) = V.dims[:left].Δ + V.dims[:right].Δ

# Implement the hashing interface (for Dict, Set)
function Base.isequal(a::Field, b::Field)
    d = (a.diagonal == b.diagonal)
    a.isKac && b.isKac && return a.r == b.r && a.s == b.s
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
    if V.degenerate
        print(io, "Degenerate $(typeof(V)), dim = ")
        show(io, V.dims[:left])
    elseif V.diagonal
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
    if V.diagonal
        if V.isKac
            if V.degenerate
                print(io, "<$(V.r), $(V.s)>")
            else
                print(io, "($(V.r), $(V.s))")
            end
        else
            P = V.dims[:left].P
            re_str = @sprintf("%.2f", real(P))
            im_str = @sprintf("%.2f", imag(P))
            if iszero(imag(P))
                print(io, "(P=$re_str)")
            elseif iszero(real(P))
                print(io, "(P=$(im_str)im)")
            else
                print(io, "(P=$(re_str)+$(im_str)im)")
            end
        end
    else
        print(io, "$(indices(V))")
    end
end
