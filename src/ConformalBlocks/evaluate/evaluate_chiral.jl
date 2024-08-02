function series_argument_from_parameter(x, s::Symbol)
    s === :crossratio && return 16*qfromx(x)
    s === :τ && return exp(im*oftype(x, π)*x)
    s === :q && return x*16
    error("the symbol you provided for the position parameter is not supported: choose
    between :crossratio, :τ, :q")
end

function evaluate_series(b::Block, x, s::Symbol, lr, der=false, reg=false)
    q = series_argument_from_parameter(x, s)
    return evalpoly(q, b._coefficients[lr])
end

function get_position(x, b::FourPointBlock)
    return crossratio(b.channel, x)
end

function get_position(x, b::OnePointBlock)
    return x
end

"""
    evaluate_chiral(b::Block, pos, s::Symbol=:crossratio)

Evaluate a block, including prefactors. If `s` is set to `:crossratio`, then pos is 
taken as the crossratio `x`, if `s` is set to :modular then `pos` is interpreted as the 
parameter `q`.
"""
function evaluate_chiral(
    c::CentralCharge,
    corr::Correlation,
    b::Block,
    pos::Number, 
    lr;
    der=false,
    reg=false
)
    x = get_position(pos, b)
    h = evaluate_series(b, x, :crossratio, lr, der, reg)
    p = blockprefactor_chiral(c, corr, b, x, lr)

    return p * h
end