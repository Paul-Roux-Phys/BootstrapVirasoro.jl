include("chiral_blocks.jl")
include("non_chiral_blocks.jl")
include("linear_combination_blocks.jl")
include("evaluate.jl")

function Block(co::CCo, chan::Symbol, d::CD=CD(), Δmax=nothing; der=false)
    Δmax === nothing && (Δmax = co.Δmax)
    CBlock(co, chan, d, Δmax, der)
end

function Block(co::CCo, d::CD=CD(), Δmax=nothing; der=false)
    Δmax === nothing && (Δmax = co.Δmax)
    CBlock(co, :s, d, Δmax, der)
end

function Block(co::NCCo, chan::Symbol, V::Field=CD(), Δmax=nothing; interchiral=false, itcr=false)
    Δmax === nothing && (Δmax = co.Δmax)
    (interchiral || itcr) ? IBlock(co, chan, V, Δmax) : NCBlock(co, chan, V, Δmax)
end

function Block(co::NCCo, V::Field=CD(), Δmax=nothing; interchiral=false, itcr=false)
    Δmax === nothing && (Δmax = co.Δmax)
    (interchiral || itcr) ? IBlock(co, :τ, V, Δmax) : NCBlock(co, :τ, V, Δmax)
end

# function Block(
#     co::Correlation{T,U},
#     chan,
#     d,
#     lr = nothing;
#     Δmax::Union{Int,Nothing} = nothing,
#     Δmax = nothing,
#     interchiral = false,
#     s_shift = 2,
#     parity = 0, # 0 for no reflection, ±1 for even/odd
# ) where {T,U}
#     # compute Δmax for this block
#     if Δmax === nothing && Δmax === nothing
#         Δmax = co.Δmax
#     elseif Δmax === nothing && Δmax !== nothing
#         Δmax = N_max(d, Δmax)
#     end
#     # don't exceed the N for which we computed the residues
#     Δmax = min(Δmax, co.Δmax)
#     if interchiral
#         @assert Δmax !== nothing "must provide Δmax for interchiral block"
#         b = InterchiralBlock(co, chan, d, Δmax, s_shift)
#     else
#         if lr === nothing
#             b = Block(co, chan, d, Δmax)
#         else
#             b = Block(co, chan, d, lr, Δmax)
#         end
#     end
#     (parity == 0 || real(d.s) <= 0 || islogarithmic(b) || isdegenerate(b)) && return b
#     b_refl = reflect(b)
#     return LinearCombinationBlock([b, b_refl], [one(T), parity*one(T)])
# end

# Evaluate blocks with b(z)
function (b::Block)(args...)
    eval(b, args...)
end
