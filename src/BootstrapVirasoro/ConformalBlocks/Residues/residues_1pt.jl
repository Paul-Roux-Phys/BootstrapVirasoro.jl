#=
# This file computes the residues of conformal blocks for one-point torus correlations,
# for use in the Zamolodchikov recursion.
=#

function Rmn_term_vanishes(r, s, D::OneDimension)
    d = D[1]
    # The term r, s in Rmn is zero if r + i = 0 and s + j = 0
    if (d.isKac && d.r + r == 0) && d.r + s == 0
        return true
    end

    return false
end

function Rmn_zero_order(m, n, D::OneDimension)
    d = D[1]
    if d.isKac && d.r%2==1 && d.s%2==1 && abs(d.r) <= 2*m-1 && abs(d.s) <= 2*n-1
        return 1
    end
    return 0
end

function Rmn_term_nonzero(r, s, d::FourDimensions)
    B = d[1].c.B
    δ1 = d[1].δ
    return δrs(r, s, B) - δ1 
end

Rmn_term_reg(r, s, d::OneDimension) = -2*d[1].P

function Rmn_term(r, s, d::OneDimension)
    Rmn_term_vanishes(r, s, d) && return Rmn_term_reg(r, s, d)
    return Rmn_term_nonzero(r, s, d)
end

function computeRmn(m::Int, n::Int, d::OneDimension)
    B = d[1].c.B
    res = prod(Rmn_term(r, s, d) for r in 1:2:2*m-1 for s in 1-2n:2:2n-1)
    return res / (2 * Dmn(m, n, B))
end