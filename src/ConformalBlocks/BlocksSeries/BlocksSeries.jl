#=
# This file defines a function H that computes the coefficients of the series expansion 
# of one-point and four-point conformal blocks.
=#


include("four_point.jl")
include("one_point.jl")

"""
Term of order N in the expansion of the zamolodchikov H-series expansion of a conformal
block.
If reg=true, for the regularised block instead.
If der=true, for the derivative of the block instead 

TODO: the if condition when reg=true is a bit clunky, could be cleaned up
"""
function series_H_N(N, B, b::Block{T}, lr, der, reg) where {T}
    (N == 0 && !der) && return one(T)
    (N == 0 &&  der) && return zero(T)
    V = b.channel_field
    P = V.P[lr]
    CN(m, n) = b._CNmn[(N, m, n)][lr]
    ind = [(m, n) for m in 1:N for n in 1:N if (N, m, n) in keys(b._CNmn)]
    isempty(ind) && return 0
    function coeff(m, n)
        if der
            -2P / (P^2 - δrs(m, n, B))^2
        elseif reg && V.isKac && V.r % 1 == V.s % 1 == 0 && V.r > 0 &&
               (lr == left && V.s > 0 || lr == right && V.s < 0) &&
               (m, n) == (V.r, abs(V.s))
            inv(4 * δrs(V.r, abs(V.s), B))
        else
            inv(P^2 - δrs(m, n, B))
        end
    end
    return sum(CN(m, n) * coeff(m, n) for (m, n) in ind)
end

function series_H(
    c::CentralCharge, b::Block, lr, der=false, reg=false
)
    return [series_H_N(N, c.B, b, lr, der, reg) for N in 0:b.Nmax]
end

