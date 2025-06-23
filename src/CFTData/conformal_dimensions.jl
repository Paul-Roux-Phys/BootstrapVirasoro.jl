"""
    ConformalDimension(c; param=value, r=0, s=0)
    
Type representing a conformal dimension, with precision given by the central charge.
The supported parameters are `Δ`, `δ`, `P`, `p`, or the Kac indices `r, s`.

# Examples

```jldoctest
julia> c = CentralCharge(β = 0.3im);

julia> d1 = ConformalDimension(c, δ=0.5); d1.Δ ≈ 3.800277777777777 + 0.0im
true

julia> d1.P ≈ 0.7071067811865476
true

julia> d2 = ConformalDimension(c, r=3//2, s=2//3)
ConformalDimension{ComplexF64} with Kac indices r = 3//2, s = 2//3
```
"""
struct ConformalDimension{T}
    c::CentralCharge{T}
    P::T
    p::T
    δ::T
    Δ::T
    r::Union{Rational,Int}
    s::Union{T,Rational,Int}
    isKac::Bool
end

# convenience alias
const CD = ConformalDimension

function Pfrom(s::Symbol, x, c::CentralCharge)
    s === :Δ && return sqrt(complex(x - (c.c-1)/24))
    s === :δ && return sqrt(complex(x))
    s === :P && return x
    s === :p && return im*x
end

function Pto(s::Symbol, x, c::CentralCharge)
    s === :Δ && return x^2 + (c.c-1)/24
    s === :δ && return x^2
    s === :P && return x
    s === :p && return -im*x
    s === :w && return -2*cos(oftype(c.β, π)*c.β*x)
end

P_rs(r, s, β::Number) = 1/2 * (β*r - s/β)
"""
        P_rs(r, s, β)
        P_rs(r, s, c)

``P_{(r, s)} = 1/2 * (β r - s / β)``.
"""
P_rs(r, s, c::CentralCharge) = P_rs(r, s, c.β)
δrs(r, s, B::Number) = -1/4 * (B*r^2 + 2*r*s + s^2/B)
δrs(r, s, c::CentralCharge) = -1/4 * (c.B*r^2 + 2*r*s + s^2/c.B)

function ConformalDimension(
    c::CentralCharge{T},
    sym::Symbol,
    P;
    r = missing,
    s = missing,
) where {T}
    β = c.β
    if (r !== missing && s !== missing)
        r % 1 == 0 ? r = Int(r) : nothing;
        s % 1 == 0 ? s = Int(s) : nothing;
        P = (r*β-s/β)/2
        isKac = true
    else
        P = Pto(:P, Pfrom(sym, P, c), c)
        isKac = false
        r = 0
        s = 2β * P
    end
    p = Pto(:p, P, c)
    δ = Pto(:δ, P, c)
    Δ = Pto(:Δ, P, c)
    ConformalDimension{T}(c, P, p, δ, Δ, r, s, isKac)
end

function ConformalDimension(
    c::CentralCharge;
    r = missing,
    s = missing,
    Δ = missing,
    δ = missing,
    P = missing,
    p = missing,
)
    (r !== missing && s !== missing) &&
        return ConformalDimension(c, :Δ, 0, r = r, s = s)
    Δ !== missing && return ConformalDimension(c, :Δ, Δ)
    δ !== missing && return ConformalDimension(c, :δ, δ)
    P !== missing && return ConformalDimension(c, :P, P)
    p !== missing && return ConformalDimension(c, :p, p)
    return ConformalDimension(c, :Δ, 0)
end

ConformalDimension() = ConformalDimension(CentralCharge())

indices(d::CD) = d.r, d.s

function Base.:+(d1::CD, d2::CD)
    ConformalDimension(d1.c, :Δ, d1.Δ + d2.Δ)
end

"""
        shift(d, i)
Shift a conformal dimension by s -> s+i or P -> P + i/β
"""
function shift(d::CD, shift, index = :s)
    c = d.c
    if index === :r
        if d.isKac
            return ConformalDimension(c, r = d.r + shift, s = d.s)
        else
            return ConformalDimension(c, P = d.P + shift / 2 * c.β)
        end
    else
        if d.isKac
            return ConformalDimension(c, r = d.r, s = d.s + shift)
        else
            return ConformalDimension(c, P = d.P - shift / 2 / c.β)
        end
    end
end

total_dimension(d::CD) = d.Δ
isdegenerate(d::CD) = d.isKac && d.r isa Int && d.s isa Int && d.r > 0 && d.s > 0

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
