function shift(V::Field, k=1)
    if V.r != 0 # the field is non-diagonal => shift s
        return Field(V.c, Kac=true, r=V.r, s=V.s+k)
    else # the field is diagonal => shift P
        P = V.P[:left]
        return Field(V.c, P=P+k/β, diagonal=true)
    end
end

# ratio C_{r, s-1}/C_{r, s+1}
function shift_C(V1, V2, r, s)
    β = V1.c.β
    r1, s1 = V1.indices
    r2, s2 = V2.indices

    res = (-1)^(2*r2 + max(2*r, 2*r1, 2*r2, r+r1+r2))
    res *= β^(-4*β^2*s)
    res *= prod(
        gamma(1//2 + abs(pm1*r1 + pm2*r2 - r) + inv(2*β^2)*(pm1*s1 + pm2*s2 - s)) /
        gamma(1//2 + abs(pm1*r1 + pm2*r2 + r) + inv(2*β^2)*(pm1*s1 + pm2*s2 + s))
        for pm1 in (-1, 1) for pm2 in (-1, 1)
    )

    return res
end

# ratio B_{r, s-1}/B_{r, s+1}
function shift_B(r, s)
    res = (-1)^(2*r)
    res *= β^(-8*inv(β^2)*s)
    res *= prod(
        gamma(a + r - s/β^2) / gamma(a + r + s/β^2)
        for a in (0, 1, inv(β^2), 1-inv(β^2))
    )
end

# ratio D_{r, s-1}/D_{r, s+1} for four-point structure constants
shift_D(Vs::FourFields, r, s) =
    shift_C(Vs[1], Vs[2], r, s)*shift_C(Vs[3], Vs[4], r, s) / shift_B(r, s)

struct BlockBulkInterchiral{T} <: Block{T}

    r::Rational
    s::Number
    blocks::Vector{BlockNonChiral{T}}

end

function evaluate(b::BlockBulkInterchiral, x)
    r0, s0 = b.r, b.s
    res = 0
    for f in b.blocks
        r, s = f.channel_field.indices
        n = Int(s-s0)
        res += prod(
            inv(shift_D(f.fields, r0, s0+i+1))
            for i in 0:2:n-2
        ) * evaluate(f, x)
    end
end