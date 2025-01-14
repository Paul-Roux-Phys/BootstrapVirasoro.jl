realify(x) = isreal(x) ? real(x) : x

const dimension_parameter_list = (:Δ, :δ, :P, :p, :w)

"""Get P from any given parameter"""
function Pfrom(s::Symbol, x, c::CentralCharge)
    s === :Δ && return realify(sqrt(complex(x - (c.c-1)/24)))
    s === :δ && return realify(sqrt(complex(x)))
    s === :P && return x
    s === :p && return im*x
end

"""Get all parameters from P"""
function Pto(s::Symbol, x, c::CentralCharge)
    s === :Δ && return x^2 + (c.c-1)/24
    s === :δ && return x^2
    s === :P && return x
    s === :p && return -im*x
    s === :w && return -2*cos(oftype(c.β, π)*c.β*x)
end

"""
    ConformalDimension{T}

Type for encoding a conformal dimension.
The supported parameters are `Δ`, `δ`, `P`, `p`, `w`, or the Kac indices `r` and `s`.

# Examples
```jldoctest
julia> c = CentralCharge(:c, 0.5);

julia> d1 = ConformalDimension(c, P = 1.2+0.1im)
ConformalDimension{ComplexF64} with
Δ = 1.4091666666666667 + 0.24im, P = 1.2 + 0.1im

julia> d2 = ConformalDimension(c, r=2, s=3//2)
ConformalDimension{ComplexF64} with Kac indices r = 2//1, s = 3//2

julia> d1.P
1.2 + 0.1im

julia> d2.Δ
-0.020833333333333332 + 0.0im

julia> d1 + d2
ConformalDimension{ComplexF64} with
Δ = 1.3883333333333332 + 0.24im, P = 1.1913489935345947 + 0.10072615216131914im

julia> d1.w
2.05874441299789 - 0.06772184182090507im
```
"""
struct ConformalDimension{T}

    c::CentralCharge{T}
    P::T
    isKac::Bool
    r::Union{Missing, Rational}
    s::Union{Missing, Rational}

end

Prs(r, s, β::Number) = 1/2 * (β*r - s/β)
Prs(r, s, c::CentralCharge) = Prs(r, s, c.β)
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

"""
    ConformalDimension(c, parameter, value; Kac=false, r=0, s=0)
    
Constructor function for the `ConformalDimension` type.

# Examples

```jldoctest
julia> c = CentralCharge(β = 0.3im)
c = 80.20666666666665 + 0.0im, β = 0.0 - 0.3im

julia> d1 = ConformalDimension(c, δ = 0.5)
ConformalDimension{ComplexF64} with
Δ = 3.800277777777777 + 0.0im, P = 0.7071067811865476

julia> d2 = ConformalDimension(c, r=3//2, s=2//3)
ConformalDimension{ComplexF64} with Kac indices r = 3//2, s = 2//3
```
"""
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
        if r % 1 == s % 1 == 0
            return Int.((r, s))
        else
            return r, s
        end
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

function Base.show(io::IO, d::ConformalDimension{T}) where {T}
    if d.isKac
        print(io, "ConformalDimension{$T} with Kac indices r = $(d.r), s = $(d.s)")
    else
        print(io, "ConformalDimension{$T} with\nΔ = $(d.Δ), P = $(d.P)")
    end
end

function Base.:+(d1::ConformalDimension, d2::ConformalDimension)
    ConformalDimension(d1.c, :Δ, d1.Δ + d2.Δ)
end

function shift(d::ConformalDimension, shift)
    c = d.c
    if d.isKac
        return ConformalDimension(c, r=d.r, s=d.s + shift)
    else
        return ConformalDimension(c, P=d.P + shift / c.β)
    end
end