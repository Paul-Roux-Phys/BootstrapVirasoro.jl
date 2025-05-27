struct BlockInterchiral{T} <: Block{T}

    blocks::Vector{Block{T}}
    shifts::Vector{T}

end

function my_mod2(s)
    s = s % 2
    return s > 1 ? s - 2 : s
end

function BlockInterchiral(
    co::Correlation{T}, chan, V, Δmax, s_shift=2
) where {T}
    @assert real(V.c.β^2) > 0 "interchiral blocks are defined for Re(β^2) > 0"
    if s_shift == 1
        error("Interchiral blocks with shifts s->s+1 not implemented")
    end
    fields = []
    shifts = []
    V0 = V
    if !isdiagonal(V) && V0.s != 1 # bring back s in (-1, 1)
        V0 = Field(V.c, r=V.r, s=my_mod2(V.s))
    end
    Vshift = V0
    D = 1
    ext_Vs = permute_fields(co.fields, chan)
    while real(Vshift.Δ[:left] + Vshift.Δ[:right]) <= real(Δmax)
        push!(fields, Vshift)
        push!(shifts, D)
        D /= shift_D(ext_Vs, Vshift)
        Vshift = shift(Vshift, s_shift)
    end
    Vshift = shift(V0, -s_shift)
    D = shift_D(ext_Vs, Vshift)
    while real(Vshift.Δ[:left] + Vshift.Δ[:right]) <= real(Δmax)
        push!(fields, Vshift)
        push!(shifts, D)
        Vshift = shift(Vshift, -s_shift)
        D *= shift_D(ext_Vs, Vshift)
    end
    blocks = [Block(co, chan, V, Δmax=Δmax) for V in fields]
    BlockInterchiral{T}(blocks, shifts)
end

""""
        shift_C(V1, V2, V3)

Ratio ``C_{123}/C_{123^++}``.
See arXiv:2411.17262 (4.13a).
"""
function shift_C123(V1, V2, V3)
    β = V1.c.β
    r1, s1 = V1.indices
    r2, s2 = V2.indices
    r3, s3 = V3.indices
    s3 += 1

    res = (-1)^(Int(2r2 + max(2r1, 2r2, 2r3, r1 + r2 + r3)))
    res *= (β+0im)^(-4inv(β)^2 * s3)
    res *= prod(
        gamma(1 // 2 + 1 // 2 * abs(pm1 * r1 + pm2 * r2 - r3) +
              inv(2 * β^2) * (pm1 * s1 + pm2 * s2 - s3)) /
        gamma(1 // 2 + 1 // 2 * abs(pm1 * r1 + pm2 * r2 + r3) +
              inv(2 * β^2) * (pm1 * s1 + pm2 * s2 + s3))
        for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

"""
        shift_C122(V1, V2)

Ratio ``C_{122}/C_{12^++2^++}``.
"""
function shift_C122(V1, V2)
    β = V1.c.β
    r1, s1 = V1.indices
    r2, s2 = V2.indices
    s2 += 1

    res = (-1)^(2r2)
    res *= β^(-8inv(β)^2 * s2)
    res *= prod(
        gamma(1 // 2 + 1 // 2 * abs(r1 + 2pm2 * r2) -
              inv(β^2) / 2 * (pm2 * s1 + 2s2 + pm1)) /
        gamma(1 // 2 + 1 // 2 * abs(r1 + 2pm2 * r2) +
              inv(β^2) / 2 * (pm2 * s1 + 2s2 + pm1))
        for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

""""
        shift_B(V)

Ratio ``B_{V}/B_{V^++}``.
See arXiv:2411.17262 (4.13b).
"""
function shift_B(V)
    β = V.c.β
    r, s = V.indices
    s += 1
    res = (-1)^(2r)
    res *= β^(-8inv(β)^2 * s)
    res *= prod(
        gamma(a + r - s / β^2) / gamma(a + r + s / β^2)
        for a in (0, 1, inv(β^2), 1 - inv(β^2))
    )
end

"""ratio ``D_{V}/D_{V^++}`` for four-point structure constants"""
function shift_D(Vs::FourFields, V)
    shift_C123(Vs[1], Vs[2], V)*shift_C123(Vs[3], Vs[4], V) / shift_B(V)
end

""" ratio ``D_{V}/D_{V^++}`` for torus one-point structure constants"""
function shift_D(Vs::OneField, V)
    shift_C122(Vs[1], V) / shift_B(V)
end

function Base.getproperty(b::BlockInterchiral, s::Symbol)
    s === :fields && return [block.channel_field for block in getfield(b, :blocks)]
    s === :indices && return [block.channel_field.indices for block in getfield(b, :blocks)]
    s in (:r, :s, :channel, :channel_field) && return getproperty(getfield(b, :blocks)[1], s)
    s === :shift && return length(b.blocks) > 1 ?
        Int(abs(getfield(b, :blocks)[2].s - getfield(b, :blocks)[1].s)) : 2
    getfield(b, s)
end

function Base.length(b::BlockInterchiral)
    length(b.fields)
end

function Base.show(io::IO, ::MIME"text/plain", b::BlockInterchiral)
    svals = [V.s for V in b.fields]
    V = b.fields[1]
    ns = real.(round.(V.s .- svals))
    shift_range = Int(minimum(ns)):b.shift:Int(maximum(ns))
    if isempty(b.fields)
        print(io, "Empty interchiral block")
    elseif isdiagonal(b.fields[1])
        println(io, "Interchiral block with channel $(b.channel) and fields")
        print(io, "V_{P = $(V.P[:left]) + n/β}, n ∈ $(shift_range)}")
    else
        println(io, "Interchiral block with channel $(b.channel) and fields")
        print(io, "{ V_{$(b.r), $(b.s) + s}, s ∈ $(shift_range) }")
    end
end

function Base.show(io::IO, b::BlockInterchiral)
    svals = [V.s for V in b.fields]
    V = b.fields[1]
    ns = real(round.(V.s .- svals))
    shift_range = Int(minimum(ns)):b.shift:Int(maximum(ns))
    if isempty(b.fields)
        print(io, "Empty interchiral block")
    elseif b.fields[1].r == 0
        print(io, "G^(s)({ ")
        print(io, "V_{P = $(V.P[:left]) + n/β}, n ∈ $(shift_range)}")
    else
        print(io, "G^(s)({ ")
        print(io, "{ V_{$(b.r), $(b.s) + s}, s ∈ $(shift_range) }")
    end
end
