"""
        LeftRight{T}
Left and right pairs of objects. Can be accessed with
`obj[:left]` and `obj[:right]`.
"""
const LeftRight{T} = Tuple{T,T} # left and right pairs of things

function Base.getindex(x::LeftRight, s::Symbol)
    s === :left && return x[1]
    s === :right && return x[2]
    error("tried to access pair element other than 1, 2, :left or :right")
end

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
    isKac::Bool
    r::Union{Missing,Rational}
    s::Union{Missing,Rational}

end

# convenience alias
const CD = ConformalDimension

"""
    Field(c, parameter=value; r, s, diagonal=false)

Type for representing a non-chiral field.

# keyword arguments:

- `r::Rational`,`s::Rational`. By convention ``V_{(r,s)}`` has left and right momenta ``(P_{(r,s)}, P_{(r,-s)})``.
- `diagonal::Bool`: set to `true` to get a diagonal field;

# Examples

```jldoctest
julia> c = CentralCharge(β = big"0.5");

julia> V = Field(c, r=0, s=1)
Diagonal Field{Complex{BigFloat}}, dim = Δ_{0, 1}

julia> V.Δ
(0.4375 + 0.0im, 0.4375 + 0.0im)

julia> V.P[:left] ≈ 1
true

julia> V.p[:right] ≈ -im
true

julia> V2 = Field(c, :P, 0.42, diagonal=true); isdiagonal(V2)
true
```
"""
struct Field{T}

    dims::LeftRight{ConformalDimension{T}}
    diagonal::Bool

end

include("CFTData/central_charges.jl")
include("CFTData/conformal_dimensions.jl")
include("CFTData/fields.jl")
