"""
    CentralCharge(param = value)

Type representing a central charge, with precision given by the type of `value`.
The supported `param` are `c`, `β`, `b`, `B`.

# Examples

```jldoctest
julia> c = CentralCharge(β=big"0.1"+big"0.2"*im); c.β ≈ big"0.1" + big"0.2"*im
true

julia> c = CentralCharge(c = 0.7);

julia> c.b ≈ -0.0 + 0.894427190999916im
true

julia> c.n ≈ 1.6180339887498953 + 0.0im
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

function CentralCharge(s::Symbol, x)
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
    β !== missing && return CentralCharge(:β, β)
    c !== missing && return CentralCharge(:c, c)
    B !== missing && return CentralCharge(:B, B)
    b !== missing && return CentralCharge(:b, b)
    n !== missing && return CentralCharge(:n, n)
    return CentralCharge(:c, 1)
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
