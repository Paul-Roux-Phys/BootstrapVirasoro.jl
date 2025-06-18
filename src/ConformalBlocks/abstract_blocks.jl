"""
        Block(co, chan, d or V; Δmax=10., interchiral=false, der=false)

Object representing a block in a channel. Can be a chiral, derivative, non-chiral, or interchiral block, depending on the types of the input correlation and channel dimension `d` or channel field `V`. Can be evaluated at a position with [`evaluate`](@ref).

# Example

```jldoctest
julia> c = CentralCharge(β = sqrt(2));

julia> V1 = Field(c, r=2, s=3//2);

julia> co = Correlation(V1, V1, V1, V1, 10);

julia> V = Field(c, r=2, s=1//2); b = Block(co, :s, V)
Non chiral factorised block for the correlation
< (2, 3//2) (2, 3//2) (2, 3//2) (2, 3//2) >
channel: s, V_{(2, 1//2)}

julia> b = Block(co[:left], :s, V.dims[:left])
Chiral block for the correlation
< Δ_{2, 3//2} Δ_{2, 3//2} Δ_{2, 3//2} Δ_{2, 3//2} >
Channel: s, Δ_{2, 1//2}
```
"""
abstract type Block{T, U} end # general conformal block. Can be interchiral, non-chiral or chiral

N_max(d::CD, Δmax) = max(0, ceil(Int, real(Δmax - d.Δ)))
N_max(V::Field, Δmax) = max(0, ceil(Int, real(Δmax - total_dimension(V))))

include("correlations.jl")
include("chiral_blocks.jl")
include("non_chiral_blocks.jl")
include("linear_combination_blocks.jl")
include("evaluate.jl")

Block(co::Corr, chan::Symbol, d::CD, Nmax::Int) = ChiralBlock(co, chan, d, Nmax)
Block(co::Corr, chan::Symbol, V::Field, Nmax::Int) = NonChiralBlock(co, chan, V, Nmax)
Block(co::Corr, chan, V::Field, lr::Symbol, Nmax; der = false) =
    ChiralBlock(co, chan, V.dims[lr], lr, Nmax, der = der)
Block(co, chan, d::CD, lr, Nmax; der = false) =
    ChiralBlock(co, chan, d, lr, Nmax, der = der)
Block(bs::Vector{<:Block}, coeffs) = LinearCombinationBlock(bs, coeffs)

function Block(
    co::Correlation{T,U},
    chan,
    d,
    lr = nothing;
    Nmax::Union{Int,Nothing} = nothing,
    Δmax = nothing,
    interchiral = false,
    s_shift = 2,
    parity=0 # 0 for no reflection, ±1 for even/odd
) where {T,U}
    # compute Nmax for this block
    if Nmax === nothing && Δmax === nothing
        Nmax = co.Nmax
    elseif Nmax === nothing && Δmax !== nothing
        Nmax = N_max(d, Δmax)
    end
    # don't exceed the N for which we computed the residues
    Nmax = min(Nmax, co.Nmax)
    if interchiral
        @assert Δmax !== nothing "must provide Δmax for interchiral block"
        b = InterchiralBlock(co, chan, d, Δmax, s_shift)
    else
        if lr === nothing
            b = Block(co, chan, d, Nmax)
        else
            b = Block(co, chan, d, lr, Nmax)
        end
    end
    (parity == 0 || d.s <= 0 || islogarithmic(b) || isdegenerate(b)) && return b
    b_refl = reflect(b)
    return LinearCombinationBlock([b, b_refl], [one(T), parity*one(T)])
end

# Evaluate blocks with b(z)
function (b::Block)(args...)
    evaluate(b, args...)
end

isdegenerate(b::Block) = isdegenerate(b.channel_field)
isdiagonal(b::Block) = isdiagonal(b.channel_field)
