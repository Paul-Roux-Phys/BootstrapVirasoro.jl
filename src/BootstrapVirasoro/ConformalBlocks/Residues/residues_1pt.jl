#=
# This file computes the residues of conformal blocks for one-point torus correlations,
# for use in the Zamolodchikov recursion.
=#

function Rmn_zero_order(m, n, V::OneField)
    v = V[1]
    if v.isKac && v.r%2==1 && v.s%2==1 && abs(v.r) <= 2*m-1 && abs(v.s) <= 2*n-1
        return 1
    end
    return 0
end

function computeRmn(m::Int, n::Int, V::OneField, lr)
    B = V[1].c.B
    δ1 = V[1].δ[lr]
    if Rmn_zero_order(m, n, V) > 0
        return 0
    else
        res = prod(prod(δrs(r, s, B) - δ1 for r in 1:2:2*m-1) for s in 1-2n:2:2n-1)
        return res/(2*Dmn(m, n, B))
    end
end

