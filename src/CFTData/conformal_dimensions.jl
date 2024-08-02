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
Type for encoding a conformal dimension, and conveniently access its values in all parametrisations
"""
struct ConformalDimension{T <: Union{AbstractFloat, Complex{Float64}, Complex{BigFloat}}}

    c::CentralCharge{T}
    P::T
    isKac::Bool
    r::Rational
    s::Rational

end

function ConformalDimension(c::CentralCharge{T}, sym::Symbol=:P, P=0; Kac=false, r=0, s=0) where {T}
    if Kac
        P = (r*c.β-s/c.β)/2
    else
        P = Pto(:P, Pfrom(sym, P, c), c)
    end
    ConformalDimension{T}(c, P, Kac, r, s)
end

function Base.getproperty(d::ConformalDimension, s::Symbol)
    c = getfield(d, :c)
    P = Pto(:P, Pfrom(:P, getfield(d, :P), c), c)
    P = realify(P)
    s in dimension_parameter_list && return Pto(s, P, c)
    return getfield(d, s)
end
