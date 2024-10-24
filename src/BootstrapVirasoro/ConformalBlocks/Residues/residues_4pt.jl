#=
# This file computes the residues of conformal blocks for four-point sphere correlations,
# for use in the Zamolodchikov recursion.
=#


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

function Rmn_term(r, s, V::FourFields, lr)
    B = V[1].c.B
    δ = [V[i].δ[lr] for i in 1:4]
    if r != 0 || s != 0
        return (((δ[2]-δ[1])^2 - 2*δrs(r, s, B)*(δ[1]+δ[2]) + δrs(r, s, B)^2)
                *((δ[3]-δ[4])^2 - 2*δrs(r, s, B)*(δ[3]+δ[4]) + δrs(r, s, B)^2))
    else
        return (δ[2]-δ[1])*(δ[3]-δ[4])
    end
end

function Rmn_term_reg(r, s, V::FourFields, lr)
    c = V[1].c
    if r != 0 || s != 0
        return 8*V[1].P[lr]*V[2].P[lr]*Field(c, Kac=true, r=r, s=s)
    else
        return 2*V[2].P[lr]
    end
end

function computeRmn(m, n, V::FourFields{T}, lr) where {T}
    if Rmn_zero_order(m, n, V) == 0
        if m == 1
            res = prod(Rmn_term(0, s, V, lr) for s in 1-n:2:0)
        else # m > 1
            res = prod(prod(Rmn_term(r, s, V, lr)
                            for s in 1-n:2:n-1) for r in 1-m:2:-1)
            if m%2 == 1 # m odd -> treat r=0 term separately
                res *= prod(Rmn_term(0, s, V, lr) for s in 1-n:2:0)
            end
        end
    else
        # if m == 1
            res = zero(T)
        # end
    end

    return res/(2*Dmn(m, n, V[1].c.B))
end

function computeRmnreg(m, n, V::FourFields{T}, lr) where {T}
    
end

