struct InterchiralBlock{T} <: Block{T}
    Δmax::Int
    chan_field::Field{T}
    fields::Vector{Field{T}}
    indices::Vector
    corr::Corr
    blocks::Vector{Block{T}}
    shifts::Vector{T}
    chan::Symbol
end
const IBlock = InterchiralBlock

struct LinearCombinationBlock{T} <: Block{T}
    chan_field::Field{T}
    blocks::Vector{Block{T}}
    coeffs::Vector{T}
    chan::Symbol
end
const LCBlock = LinearCombinationBlock

function my_mod2(s)
    s = s % 2
    return s > 1 ? s - 2 : s
end

function LinearCombinationBlock(bs::Vector{<:Block{T}}, coeffs) where {T}
    LinearCombinationBlock{T}(bs[1].chan_field, bs, coeffs, bs[1].chan)
end

function Base.:+(b1::LCBlock, b2::LCBlock)
    LCBlock(vcat(b1.blocks, b2.blocks), vcat(b1.coeffs, b2.coeffs))
end

function Base.:+(b1::LCBlock, b2::Block)
    LCBlock(vcat(b1.blocks, b2), vcat(b1.coeffs, 1))
end

function Base.:+(b1::Block, b2::LCBlock)
    LCBlock(vcat(b1, b2.blocks), vcat(1, b2.coeffs))
end

function Base.:+(b1::Block, b2::Block)
    LCBlock(vcat(b1, b2), vcat(1, 1))
end

function Base.:-(b1::Block, b2::Block)
    LCBlock(vcat(b1, b2), vcat(1, -1))
end

function Base.:*(a::Number, b::Block)
    LCBlock([b], [a])
end

function InterchiralBlock(co::Co{T}, chan, V, Δmax) where {T}
    @assert real(V.c.β^2) > 0 "interchiral blocks are implemented for Re(β^2) > 0"
    fields = []
    shifts = []
    extfields = getfields(co, chan)
    V0 = V
    if !V.diagonal
        V0 = Field(V.c, r = V.r, s = my_mod2(V.s))
    end
    Vshift = V0
    D = 1
    while real(total_dimension(Vshift)) <= real(Δmax) && spin(Vshift) < real(Δmax)
        push!(fields, Vshift)
        push!(shifts, D)
        D /= shift_D(extfields, Vshift)
        Vshift = shift(Vshift, 2)
    end
    if !(V.r isa Int && V.s isa Int && V.r > 0 && V.s > 0)
        Vshift = shift(V0, -2)
        D = shift_D(extfields, Vshift)
        while real(total_dimension(Vshift)) <= real(Δmax) && spin(Vshift) < real(Δmax)
            push!(fields, Vshift)
            push!(shifts, D)
            Vshift = shift(Vshift, -2)
            D *= shift_D(extfields, Vshift)
        end
    end
    blocks = [Block(co, chan, V, Δmax) for V in fields]
    ind = [indices(b.chan_field) for b in blocks]

    InterchiralBlock{T}(Δmax, V, fields, ind, co, blocks, shifts, chan)
end

function islogarithmic(b::IBlock)
    isempty(b.blocks) && return false
    islogarithmic(b.blocks[1])
end

reflect(b::IBlock) = IBlock(b.corr, b.chan, reflect(b.chan_field), b.Δmax)

""""
        shift_C(V1, V2, V3)

Ratio ``C_{123}/C_{123^++}``.
See arXiv:2411.17262 (4.13a).
"""
function shift_C123(V1, V2, V3)
    β = V1.c.β
    r1, s1 = get_indices(V1)
    r2, s2 = get_indices(V2)
    r3, s3 = get_indices(shift(V3, 1))

    res = (-1)^(Int(2r2 + max(2r1, 2r2, 2r3, r1 + r2 + r3)))
    res *= (β+0im)^(-4/β^2 * s3)
    res *= prod(
        gamma(
            1 // 2 +
            1 // 2 * abs(pm1 * r1 + pm2 * r2 - r3) +
            1/(2 * β^2) * (pm1 * s1 + pm2 * s2 - s3),
        ) / gamma(
            1 // 2 +
            1 // 2 * abs(pm1 * r1 + pm2 * r2 + r3) +
            1/(2 * β^2) * (pm1 * s1 + pm2 * s2 + s3),
        ) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

"""
        shift_C122(V1, V2)

Ratio ``C_{122}/C_{12^++2^++}``.
"""
function shift_C122(V1, V2)
    β = V1.c.β
    r1, s1 = get_indices(V1)
    r2, s2 = get_indices(shift(V2, 1))

    res = (-1)^(Int(2r2))
    res *= β^(-8/β^2 * s2)
    res *= prod(
        gamma(
            1 // 2 + 1 // 2 * abs(r1 + 2pm2 * r2) -
            1/(β^2 * 2) * (pm2 * s1 + 2s2 + pm1),
        ) / gamma(
            1 // 2 +
            1 // 2 * abs(r1 + 2pm2 * r2) +
            1/(β^2 * 2) * (pm2 * s1 + 2s2 + pm1),
        ) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

""""
        shift_B(V)

Ratio ``B_{V}/B_{V^++}``.
See arXiv:2411.17262 (4.13b).
"""
function shift_B(V)
    β = V.c.β
    r, s = get_indices(shift(V, 1))
    res = (-1)^(2r)
    res *= β^(-8inv(β)^2 * s)
    res *= prod(
        gamma(a + r - s / β^2) / gamma(a + r + s / β^2) for
        a in (0, 1, inv(β^2), 1 - inv(β^2))
    )
end

"""ratio ``D_{V}/D_{V^++}`` for four-point structure constants"""
function shift_D(Vs::NTuple{4, Field}, V)
    shift_C123(Vs[1], Vs[2], V)*shift_C123(Vs[3], Vs[4], V) / shift_B(V)
end

""" ratio ``D_{V}/D_{V^++}`` for torus one-point structure constants"""
function shift_D(Vs::NTuple{1, Field}, V)
    shift_C122(Vs[1], V) / shift_B(V)
end

# function Base.getproperty(b::InterchiralBlock, s::Symbol)
#     s === :indices &&
#         return [indices(block.channel_field) for block in getfield(b, :blocks)]
#     s in (:r, :s, :channel, :channel_field, :corr, :correlation) &&
#         return getproperty(getfield(b, :blocks)[1], s)
#     s === :shift && return length(b.blocks) > 1 ?
#            Int(abs(getfield(b, :blocks)[2].s - getfield(b, :blocks)[1].s)) : 2
#     getfield(b, s)
# end

# function Base.getproperty(b::LCBlock, s::Symbol)
#     s === :corr && return getfield(b, :blocks)[1].correlation
#     getfield(b, s)
# end

function Base.length(b::IBlock)
    length(b.fields)
end

# function Base.getproperty(b::IBlock, s::Symbol)
#     s in (:r, :s) && return getproperty(getfield(b, :blocks)[1], s)
#     getfield(b, s)
# end

function Base.show(io::IO, ::MIME"text/plain", b::IBlock)
    svals = [V.s for V in b.fields]
    V = b.fields[1]
    ns = real.(round.(-V.s .+ svals)) ./ 2
    shift_range = Int(minimum(ns)):1:Int(maximum(ns))
    if isempty(b.fields)
        print(io, "Empty interchiral block")
    elseif !b.fields[1].diagonal || b.fields[1].degenerate
        println(io, "Interchiral block with channel $(b.chan) and fields")
        print(io, "{ V_{$(b.chan_field.r), $(b.chan_field.s) + 2s}, s ∈ $(shift_range) }")
    else
        println(io, "Interchiral block with channel $(b.chan) and fields")
        print(io, "V_{P = $(V.dims.left.P) + n/β}, n ∈ $(shift_range)}")
    end
end

function Base.show(io::IO, b::IBlock)
    svals = [V.s for V in b.fields]
    V = b.fields[1]
    ns = real(round.(-V.s .+ svals)) ./ 2
    shift_range = Int(minimum(ns)):1:Int(maximum(ns))
    if isempty(b.fields)
        print(io, "Empty interchiral block")
    elseif !b.fields[1].diagonal || b.fields[1].degenerate 
        print(io, "G^(s)({ ")
        print(io, "V_{$(b.chan_field.r), $(b.chan_field.s) + 2s}, s ∈ $(shift_range) })")
    else
        print(io, "G^(s)({ ")
        print(io, "V_{P = $(V.dims[:left].P) + n/β}, n ∈ $(shift_range)})")
    end
end

function Base.show(io::IO, b::LCBlock)
    coeffs = [
        if c ≈ round(Int, c)
            round(Int, c)
        else
            c
        end for c in b.coeffs
    ]
    print(io, "($(coeffs[1])) * $(b.blocks[1])")
    for (i, bl) in enumerate(b.blocks[2:end])
        print(io, " + ($(coeffs[i+1])) * $bl")
    end
end
