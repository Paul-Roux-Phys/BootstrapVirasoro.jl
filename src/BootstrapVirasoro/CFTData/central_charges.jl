"""
    CentralCharge{T}

Type representing a central charge.
T is expected to be a real or complex number, of standard or arbitrary precision
"""
struct CentralCharge{T}
    β::T
end

"""Get B from given parameter"""
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
end

"""Get asked parameter from B"""
function Bto(s::Symbol, x)
    rx = sqrt(complex(x))
    s === :β && return -im*rx
    s === :c && return 13+6*x+6/x
    s === :b && return -rx
    s === :B && return x
end


function Base.getproperty(c::CentralCharge, s::Symbol)
    β = Bto(:β, Bfrom(:β, getfield(c, :β)))
    s === :β && return β
    s === :c && return 13 - 6*β^2 - 6/β^2
    s === :B && return -β^2
    s === :b && return -im*β
    s === :n && return -2*cos(oftype(β, π)*β^2)
    error("$s is not a supported parametrisation of the central charge")
end

"""
    CentralCharge(parameter, value)

Constructor function for the CentralCharge type.

Given one of the four parameters `c`, `b`, `β`, `B` and its value,
creates an object CentralCharge{T} where T is real if `β` is real.

# Example
```julia-repl
julia> setprecision(BigFloat, 20, base=10)
julia> CentralCharge(big"1.2")
c = 0.1933333333333333332741, β = 1.200000000000000000003

```
"""
function CentralCharge(s::Symbol, x)
    β = Bto(:β, Bfrom(s, x))
    CentralCharge(β)
end

function CentralCharge(; β=missing, c=missing, B=missing, n=missing, b=missing)
    β !== missing && return CentralCharge(:β, β)
    c !== missing && return CentralCharge(:c, c)
    B !== missing && return CentralCharge(:B, B)
    n !== missing && return CentralCharge(:n, n)
    b !== missing && return CentralCharge(:b, b)
    return CentralCharge(:c, 1)
end

"""Display an object of type CentralCharge"""
function Base.show(io::IO, c::CentralCharge)
    println("c = $(c.c), β = $(c.β)")
end

function Base.isequal(c1::CentralCharge, c2::CentralCharge)
    return c1.β == c2.β
end

function Base.hash(c::CentralCharge, h::UInt)
    return hash(c.β, h)
end