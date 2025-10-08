function Cref(V₁, V₂, V₃, DG)
    β = V₁.c.β
    r₁, s₁ = get_indices(V₁)
    r₂, s₂ = get_indices(V₂)
    r₃, s₃ = get_indices(V₃)

    return prod(
        1 / DG(
            (β + 1 / β) / 2 +
            β / 2 * abs(pm₁ * r₁ + pm₂ * r₂ + pm₃ * r₃) +
            1 / 2 / β * (pm₁ * s₁ + pm₂ * s₂ + pm₃ * s₃),
        ) for pm₁ in (-1, 1) for pm₂ in (-1, 1) for pm₃ in (-1, 1)
    )
end

Cref(V₁, V₂, V₃) = Cref(V₁, V₂, V₃, DoubleGamma(V₁.c.β))

function Bref(DG, c, r, s, reg = 1/big(10^15))
    β = c.β
    if r%1 == 0 && s % 1 == 0
        s += reg
    end
    π = convert(typeof(c).parameters[1], Base.π) # π in the correct precision
    return (-1)^(round(Int, r*s)) / 2 / sin(π * (r%1 + s)) / sin(π * (r + s/β^2)) /
           prod(
        DG(β + pm1 * β * r + pm2 * s / β) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

function Bref(DG, c, P)
    β = c.β
    prod(1/DG(β^pm1 + pm2 * 2P) for pm1 in (-1, 1) for pm2 in (-1, 1))
end

function Bref(V::Field, DG, reg = 1/big"10"^15)
    c = V.c
    if V.diagonal
        return Bref(DG, c, V.dims[:left].P)
    else
        return Bref(DG, c, V.r, V.s, reg)
    end
end

function compute_reference(b::Block{T,U}, DG) where {T,U<:FourPoints}
    chan = b.channel
    V₁, V₂, V₃, V₄ = permute_fields(b.correlation.fields, chan)
    V = b.channel_field
    return Cref(V₁, V₂, V, DG) * Cref(V₃, V₄, V, DG) / Bref(V, DG)
end

function compute_reference(b::Block{T,U}, DG) where {T,U<:OnePoint}
    V₁ = b.correlation.fields[1]
    V = b.channel_field
    return Cref(V₁, V, V, DG) / Bref(V, DG)
end
