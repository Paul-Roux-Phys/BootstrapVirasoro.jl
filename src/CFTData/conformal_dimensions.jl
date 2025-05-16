realify(x) = isreal(x) ? real(x) : x

const dimension_parameter_list = (:Δ, :δ, :P, :p, :w)

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

"""
    ConformalDimension(c; param=value, r=0, s=0)
    
Type representing a conformal dimension, with precision given by the central charge.
The supported parameters are `Δ`, `δ`, `P`, `p`, or the Kac indices `r, s`.

# Examples

```jldoctest
julia> c = CentralCharge(β = 0.3im);

julia> d1 = ConformalDimension(c, δ=0.5)
ConformalDimension{ComplexF64} with
Δ = 3.800277777777777 + 0.0im, P = 0.7071067811865476

julia> d2 = ConformalDimension(c, r=3//2, s=2//3)
ConformalDimension{ComplexF64} with Kac indices r = 3//2, s = 2//3
```
"""
struct ConformalDimension{T}

    c::CentralCharge{T}
    P::T
    isKac::Bool
    r::Union{Missing, Rational}
    s::Union{Missing, Rational}

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
    r=missing, s=missing
) where {T}
    if (r !== missing && s !== missing)
        P = (r*c.β-s/c.β)/2
        isKac = true
    else
        P = Pto(:P, Pfrom(sym, P, c), c)
        isKac = false
    end
    ConformalDimension{T}(c, P, isKac, r, s)
end

function ConformalDimension(
    c::CentralCharge;
    r=missing, s=missing,
    Δ=missing, δ=missing, P=missing, p=missing
)
    (r !== missing && s !== missing) &&
        return ConformalDimension(c, :Δ, 0, r=r, s=s)
    Δ !== missing && return ConformalDimension(c, :Δ, Δ)
    δ !== missing && return ConformalDimension(c, :δ, δ)
    P !== missing && return ConformalDimension(c, :P, P)
    p !== missing && return ConformalDimension(c, :p, p)
    return ConformalDimension(c, :Δ, 0)
end

ConformalDimension() = ConformalDimension(CentralCharge())

function get_indices(d::ConformalDimension)
    if d.isKac
        r = getfield(d, :r)
        s = getfield(d, :s)
        if r % 1 == 0
            r = Int(r)
        end
        if s % 1 == 0
            s = Int(s)
        end
        return r, s
    else
        return 0, 2 * d.c.β * d.P
    end
end

function Base.getproperty(d::ConformalDimension, s::Symbol)
    c = getfield(d, :c)
    P = Pto(:P, Pfrom(:P, getfield(d, :P), c), c)
    P = realify(P)
    s in dimension_parameter_list && return Pto(s, P, c)
    s === :indices && return get_indices(d)
    s === :r && return get_indices(d)[1]
    s === :s && return get_indices(d)[2]
    return getfield(d, s)
end

function Base.show(io::IO, ::MIME"text/plain", d::ConformalDimension{T}) where {T}
    if d.isKac
        print(io, "ConformalDimension{$T} with Kac indices r = $(d.r), s = $(d.s)")
    else
        print(io, "ConformalDimension{$T} with\nΔ = $(d.Δ), P = $(d.P)")
    end
end

function Base.show(io::IO, d::ConformalDimension{T}) where {T}
    if d.isKac
        print(io, "Δ_{$(d.r), $(d.s)}")
    else
        print(io, "Δ_{P=$(d.P)}")
    end
end

function Base.:+(d1::ConformalDimension, d2::ConformalDimension)
    ConformalDimension(d1.c, :Δ, d1.Δ + d2.Δ)
end

function shift(d::ConformalDimension, shift, index=:s)
    c = d.c
    if index === :r
        if d.isKac
            return ConformalDimension(c, r=d.r + shift, s=d.s)
        else
            return ConformalDimension(c, P=d.P + shift / 2 * c.β)
        end
    else
        if d.isKac
            return ConformalDimension(c, r=d.r, s=d.s + shift)
        else
            return ConformalDimension(c, P=d.P + shift / 2 / c.β)
        end
    end
end

total_dimension(d::ConformalDimension) = d.Δ

function Base.hash(d::ConformalDimension, h::UInt)
    return hash((d.c, d.P, d.isKac), h)
end
