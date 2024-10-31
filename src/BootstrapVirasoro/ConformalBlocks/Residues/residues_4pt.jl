#=
# This file computes the residues of conformal blocks for four-point sphere correlations,
# for use in the Zamolodchikov recursion.
=#

function Rmn_term_vanishes(r, s, i, j, V::FourFields)
    !(V[j].isKac && V[j].isKac) && return false

    rs=[V[i].r for i in 1:4]
    ss=[V[i].s for i in 1:4]


    #= The term r, s in Rmn is zero if r1 \pm r2 + r or r3 \pm r4 + r is 0, and
    s1 \pm s2 + s or s3 \pm s4 + s is 0.
    =#
    for pm in (-1, 1)
        if (V[i].isKac && V[j].isKac
            && (rs[i] + pm * rs[j] + r == 0)
            && (ss[i] + pm * ss[j] + s == 0))
            
            return true
        end
    end

    return false
end


function Rmn_zero_order(m, n, V::FourFields)
    order = 0

    if !((V[1].isKac && V[2].isKac) || (V[3].isKac && V[4].isKac))
        return 0
    end

    r=[V[i].r for i in 1:4]
    s=[V[i].s for i in 1:4]

    #= Rmn is zero if r1 \pm r2 or r3 \pm r4 is an integer in 1-m:2:m-1, and
    s1 \pm s2 or s3 \pm s4 is an integer in 1-n:2:n-1.
    equivalently, if (|r1 \pm r2| <= m-1 and r1-r2 - (m-1) % 2 == 0)
    and (|s1 \pm s2| <= n-1 and s1-s2 - (n-1) % 2 == 0)
    =#
    for pm in (-1, 1), (i, j) in ((1, 2), (3, 4))
        if (V[i].isKac && V[j].isKac
            && (abs(r[i] + pm * r[j]) <= m - 1 && (r[i] + pm * r[j] - (m - 1)) % 2 == 0)
            && (abs(s[i] + pm * s[j]) <= n - 1 && (s[i] + pm * s[j] - (n - 1)) % 2 == 0))
            
            order += 1

        end
    end

    return order
end

function Rmn_term_nonzero(r, s, i, j, V, lr)
    B = V[1].c.B
    δ = [V[i].δ[lr] for i in 1:4]
    (r != 0 || s != 0) && return (δ[j] - δ[i])^2 - 2 *
                                 δrs(r, s, B) * (δ[i] + δ[j]) + δrs(r, s, B)^2
    return (δ[j] - δ[i]) * (-1)^(j/2)
end

function Rmn_term_reg(r, s, i, j, V::FourFields, lr)
    (r != 0 || s != 0) && return 8 * V[i].P[lr] * V[j].P[lr] * Prs(r, s, V[1].c.β)
    return 2 * V[j].P[lr]
end

function Rmn_term(r, s, i, j, V::FourFields, lr)
    Rmn_term_vanishes(r, s, i, j, V) && return Rmn_term_reg(r, s, i, j, V, lr)
    return Rmn_term_nonzero(r, s, i, j, V, lr)
end

Rmn_term(r, s, V::FourFields, lr) = prod(
    Rmn_term(r, s, i, j, V, lr) for (i, j) in ((1, 2), (3, 4))
)

function computeRmn(m, n, V::FourFields{T}, lr) where {T}
    if m == 1
        res = prod(Rmn_term(0, s, V, lr) for s in 1-n:2:0)
    else # m > 1
        res = prod(prod(Rmn_term(r, s, V, lr)
                        for s in 1-n:2:n-1) for r in 1-m:2:-1)
        if m % 2 == 1 # m odd -> treat r=0 term separately
            res *= prod(Rmn_term(0, s, V, lr) for s in 1-n:2:0)
        end
    end

    return res/(2*Dmn(m, n, V[1].c.B)) 
end