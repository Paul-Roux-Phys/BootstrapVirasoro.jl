"""
    CentralCharge{T}

Type representing a central charge.
T is expected to be a real or complex number, of standard or arbitrary precision.
The supported parameters are `c`, `β`, `b`, `B`.

# Examples

```jldoctest
julia> c = CentralCharge(c = 0.7)
c = 0.7000000000000011 + 0.0im, β = -0.894427190999916 - 0.0im

julia> c.b
-0.0 + 0.894427190999916im

julia> c.n
1.6180339887498953 + 0.0im
```
"""
struct CentralCharge{T}
    β::T
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

function CentralCharge(s::Symbol, x)
    β = Bto(:β, Bfrom(s, x))
    CentralCharge(β)
end

"""
    CentralCharge(parameter = value)

Constructor function for the CentralCharge type.

Given one of the four parameters `c`, `b`, `β`, `B` and its value,
creates an object CentralCharge{T}.

# Example

```jldoctest
julia> setprecision(BigFloat, 20, base=10);

julia> CentralCharge(B = 0.5)
c = 27.999999999999996 + 0.0im, β = 0.0 - 0.7071067811865476im

julia> CentralCharge(β = 0.7)
c = -2.184897959183676 + 0.0im, β = -0.7 - 0.0im

julia> CentralCharge(c = big"0.1" + 0.2im)
c = 0.09999999999999999991326 + 0.2000000000000000111158im, β = -0.8237591041762989640376 - 0.01729590504934815486866im
```
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
