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
If der=true, for the derivative of the block instead.
"""
function series_H_N(N, B, b::Block{T}, lr, der, reg) where {T}
    if N == 0 && !der
        return one(T)
    elseif N == 0 && der
        return zero(T)
    else
        P = b.channel_field.P[lr]
        CN(m, n) = b._CNmn[(N, m, n)][lr]
        ind = [(m, n) for m in 1:N for n in 1:N if (N, m, n) in keys(b._CNmn)]
        isempty(ind) && return 0
        if der
            return sum(2 * P * CN(m, n) / (P^2 - δrs(m, n, B))^2 for (m, n) in ind)
        else
            return sum(CN(m, n) / (P^2 - δrs(m, n, B)) for (m, n) in ind)
        end
    end
end

function series_H(
    c::CentralCharge, b::Block, lr, der=false, reg=false
)
    return [series_H_N(N, c.B, b, lr, der, reg) for N in 0:b.Nmax]
end

