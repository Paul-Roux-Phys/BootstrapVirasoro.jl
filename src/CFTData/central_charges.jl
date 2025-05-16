"""
    CentralCharge(param = value)

Type representing a central charge, with precision given by the type of `value`.
The supported `param` are `c`, `β`, `b`, `B`.

# Examples

```jldoctest
julia> c = CentralCharge(β=big"0.1"+big"0.2"*im); c.c ≈ 85.18 && c.β ≈ 0.1 + 0.2im
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
end

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
end


function Base.getproperty(c::CentralCharge, s::Symbol)
    β = getfield(c, :β)
    s in (:β, :c, :B, :b) && return getfield(c, s)
    s === :n && return -2*cos(π*β^2)
    error("$s is not a supported parametrisation of the central charge")
end

function CentralCharge(s::Symbol, x)
    B = complex(Bfrom(s, x))
    β = Bto(:β, B)
    b = Bto(:b, B)
    c = Bto(:c, B)
    CentralCharge(β, B, b, c)
end

"""
"""
function CentralCharge(; β=missing, c=missing, B=missing, b=missing)
    β !== missing && return CentralCharge(:β, β)
    c !== missing && return CentralCharge(:c, c)
    B !== missing && return CentralCharge(:B, B)
    b !== missing && return CentralCharge(:b, b)
    return CentralCharge(:c, 1)
end

function Base.show(io::IO, ::MIME"text/plain", c::CentralCharge)
    print(io, "c = $(c.c), β = $(c.β)")
end

function Base.show(io::IO, c::CentralCharge)
    print(io, "c = $(c.c)")
end

function Base.isequal(c1::CentralCharge, c2::CentralCharge)
    return c1.β == c2.β
end

function Base.hash(c::CentralCharge, h::UInt)
    return hash(c.β, h)
end
