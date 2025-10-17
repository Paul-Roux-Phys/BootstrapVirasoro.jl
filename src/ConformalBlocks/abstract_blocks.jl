include("chiral_blocks.jl")
include("non_chiral_blocks.jl")
include("linear_combination_blocks.jl")
include("evaluate.jl")

function Block(co::CCo, chan::Symbol, d::CD = CD(), Δmax = nothing; der = false)
    Δmax === nothing && (Δmax = co.Δmax)
    CBlock(co, chan, d, Δmax, der)
end

function Block(co::CCo, d::CD = CD(), Δmax = nothing; der = false)
    Δmax === nothing && (Δmax = co.Δmax)
    CBlock(co, :s, d, Δmax, der)
end

function Block(
    co::NCCo,
    chan,
    V,
    Δmax = missing;
    interchiral = false,
    itcr = false,
    parity = 0,
)
    typeof(co) <: Correlation4 &&
        chan === missing &&
        error("A four-point block must have channel :s, :t, :u. Define it with
              Block(co, chan, V, Δmax)")
    Δmax === missing && (Δmax = co.Δmax)
    b = (interchiral || itcr) ? IBlock : NonChiralBlock
    if parity == 0 || V.diagonal || V.s == 1 || V.s == 0
        return b(co, chan, V, Δmax)
    else
        return b(co, chan, V, Δmax) + parity * b(co, chan, reflect(V), Δmax)
    end
end

function Block(co::NCCo, V::Field, Δmax=missing; interchiral=false, itcr=false, parity=0)
    Block(co, :τ, V, Δmax, interchiral=interchiral, itcr=itcr, parity=parity)
end

# Evaluate blocks with b(z)
function (b::Block)(args...)
    eval(b, args...)
end
