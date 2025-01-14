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

struct BlockInterchiralNonChiral{T} <: Block{T}

    blocks::Vector{BlockNonChiral{T}}

end

function BlockInterchiralNonChiral(
    c::CorrelationNonChiral{T}, chan, V::Field, Δmax
) where {T}
    fields = []
    k = 0
    while minimum(abs(shift(V, k).dims[lr].Δ) for lr in (:left, :right)) <= abs(Δmax.Δ) ||
          minimum(abs(shift(V, -k).dims[lr].Δ) for lr in (:left, :right)) <= abs(Δmax.Δ)
        push!(fields, shift(V, k))
        push!(fields, shift(V, -k))
        k += 2
    end
    blocks = [
        Block(c, chan, V, Δmax) 
        for V in fields if abs(V.dims[:left].Δ + V.dims[:right].Δ) <= abs(Δmax.Δ)
    ]
    BlockInterchiralNonChiral{T}(blocks)
end

function Base.getproperty(b::BlockInterchiralNonChiral, s::Symbol)
    s === :fields && return [block.channel_field for block in b.blocks]
    getfield(b, s)
end

function evaluate(b::BlockInterchiralNonChiral, x)
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

function Base.show(io::IO, b::BlockInterchiralNonChiral)
    println(io, "Non chiral interchiral block with channel fields")
    for V in b.fields
        println(io, V)
    end
end