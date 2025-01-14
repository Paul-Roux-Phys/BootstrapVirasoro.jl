#=
# This file computes the residues of conformal blocks for four-point sphere correlations,
# for use in the Zamolodchikov recursion.
=#

function Rmn_term_vanishes(r, s, i, j, d::FourDimensions)
    !(d[j].isKac && d[j].isKac) && return false

    rs=[d[i].r for i in 1:4]
    ss=[d[i].s for i in 1:4]


    #= The term r, s in Rmn is zero if r1 \pm r2 + r or r3 \pm r4 + r is 0, and
    s1 \pm s2 + s or s3 \pm s4 + s is 0.
    =#
    for pm in (-1, 1), pm2 in (-1, 1)
        if (d[i].isKac && d[j].isKac
            && (rs[i] + pm * rs[j] + pm2 * r == 0)
            && (ss[i] + pm * ss[j] + pm2 * s == 0))
            
            return true
        end
    end

    return false
end

function reg_signs(r, s, i, j, d::FourDimensions)
    rs=[d[i].r for i in 1:4]
    ss=[d[i].s for i in 1:4]

    for pm in (-1, 1), pm2 in (-1, 1)
        if (rs[i] + pm * rs[j] + pm2 * r == 0) &&
           (ss[i] + pm * ss[j] + pm2 * s == 0)
            return pm, pm2
        end
    end
end

function Rmn_zero_order(m, n, d::FourDimensions)
    order = 0

    if !((d[1].isKac && d[2].isKac) || (d[3].isKac && d[4].isKac))
        return 0
    end

    r=[d[i].r for i in 1:4]
    s=[d[i].s for i in 1:4]

    #= Rmn is zero if r1 \pm r2 or r3 \pm r4 is an integer in 1-m:2:m-1, and
    s1 \pm s2 or s3 \pm s4 is an integer in 1-n:2:n-1.
    equivalently, if (|r1 \pm r2| <= m-1 and r1-r2 - (m-1) % 2 == 0)
    and (|s1 \pm s2| <= n-1 and s1-s2 - (n-1) % 2 == 0)
    =#
    for pm in (-1, 1), (i, j) in ((1, 2), (3, 4))
        if (d[i].isKac && d[j].isKac
            && (abs(r[i] + pm * r[j]) <= m - 1 && (r[i] + pm * r[j] - (m - 1)) % 2 == 0)
            && (abs(s[i] + pm * s[j]) <= n - 1 && (s[i] + pm * s[j] - (n - 1)) % 2 == 0))
            
            order += 1
        end
    end

    return order
end

function Rmn_term_nonzero(r, s, i, j, d::FourDimensions)
    B = d[1].c.B
    δ = [d[i].δ for i in 1:4]
    (r != 0 || s != 0) && return (δ[j] - δ[i])^2 - 2 *
                                 δrs(r, s, B) * (δ[i] + δ[j]) + δrs(r, s, B)^2
    return (δ[j] - δ[i]) * (-1)^(j/2)
end

function Rmn_term_reg(r, s, i, j, d::FourDimensions)
    (r != 0 || s != 0) && begin
        signs = reg_signs(r, s, i, j, d)
        r0, s0 = -signs[2] .* d[i].indices .- signs[2] * signs[1] .* d[j].indices
        P = ConformalDimension(d[1].c, r=r0, s=s0).P
        return 8 * signs[1] * signs[2] * d[i].P * d[j].P * P
    end
    return 2 * d[j].P
end

function Rmn_term(r, s, i, j, d::FourDimensions)
    Rmn_term_vanishes(r, s, i, j, d) && return Rmn_term_reg(r, s, i, j, d)
    return Rmn_term_nonzero(r, s, i, j, d)
end

Rmn_term(r, s, d::FourDimensions) = prod(
    Rmn_term(r, s, i, j, d) for (i, j) in ((1, 2), (3, 4))
)

"""
    computeRmn(m, n, d::FourDimensions)

Compute the ``Δ``-residue ``R_{m, n}`` of four-point blocks with external dimensions `d`,
or regularisation thereof, so that ratios of vanishing residues are correct.
"""
function computeRmn(m, n, d::FourDimensions{T}) where {T}
    if m == 1
        res = prod(Rmn_term(0, s, d) for s in 1-n:2:0)
    else # m > 1
        res = prod(prod(Rmn_term(r, s, d)
                        for s in 1-n:2:n-1) for r in 1-m:2:-1)
        if m % 2 == 1 # m odd -> treat r=0 term separately
            res *= prod(Rmn_term(0, s, d) for s in 1-n:2:0)
        end
    end

    return res/(2*Dmn(m, n, d[1].c.B)) 
end