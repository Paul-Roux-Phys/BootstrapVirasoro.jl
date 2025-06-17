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
    if (r !== missing && s !== missing)
        P = (r*c.β-s/c.β)/2
        isKac = true
    else
        P = Pto(:P, Pfrom(sym, P, c), c)
        isKac = false
    end
    p = Pto(:p, P, c)
    δ = Pto(:δ, P, c)
    Δ = Pto(:Δ, P, c)
    ConformalDimension{T}(c, P, p, δ, Δ, isKac, r, s)
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
    (r !== missing && s !== missing) && return ConformalDimension(c, :Δ, 0, r = r, s = s)
    Δ !== missing && return ConformalDimension(c, :Δ, Δ)
    δ !== missing && return ConformalDimension(c, :δ, δ)
    P !== missing && return ConformalDimension(c, :P, P)
    p !== missing && return ConformalDimension(c, :p, p)
    return ConformalDimension(c, :Δ, 0)
end

ConformalDimension() = ConformalDimension(CentralCharge())

function get_indices(d::CD)
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

function isdegenerate(d::CD)
    return d.isKac && d.r%1 == d.s%1 == 0 && d.r > 0 && d.s > 0
end

function Base.getproperty(d::CD, s::Symbol)
    s === :indices && return get_indices(d)
    s === :r && return get_indices(d)[1]
    s === :s && return get_indices(d)[2]
    return getfield(d, s)
end

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
