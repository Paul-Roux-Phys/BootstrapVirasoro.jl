struct ChannelSpectrum{T}
    corr::Corr
    blocks::Dict{Field{T},Block{T}}
    Δmax::Int
    chan::Symbol
end
const ChanSpec = ChannelSpectrum

fields(s::ChanSpec) = keys(s.blocks)

function _add_one!(s::ChanSpec, b)
    if !(b.chan_field in fields(s)) &&
       real(total_dimension(b.chan_field)) < real(s.Δmax)
        (s.blocks)[b.chan_field] = b
        return b
    end
end

_add!(s::ChanSpec, b::Block) = _add_one!(s, b)
_add!(s::ChanSpec, blocks::Vector) = foreach(b -> _add_one!(s, b), blocks)
add!(s, bs) = _add!(s, bs)

function add(s, blocks)
    s2 = deepcopy(s)
    add!(s2, blocks)
    return s2
end

function hasdiagonals(s::ChanSpec)
    for V in keys(s.blocks)
        V.diagonal && return true
    end
    return false
end

function ChannelSpectrum(co::Co{T}, chan, Vs::Vector{<:Field}, f::Function) where {T}
    s = ChannelSpectrum{T}(co, Dict{Field{T},Block{T}}(), co.Δmax, chan)
    for V in Vs
        if real(total_dimension(V)) < s.Δmax
            add!(s, f(V))
        end
    end
    return s
end

ChannelSpectrum(co::Co, Vs::Vector{<:Field}, f::Function) =
    ChannelSpectrum(co, :τ, Vs, f)

function _remove_one!(s::ChanSpec, V)
    delete!(s.blocks, V)
end

_remove!(s::ChanSpec, V) = _remove_one!(s, V)
_remove!(s::ChanSpec, fields::Vector) = foreach(V -> _remove_one!(s, V), fields)
remove!(s, V) = _remove!(s, V)
remove(s, Vs) = remove!(deepcopy(s), Vs)

function Base.filter(f, s::ChanSpec{T}) where {T}
    return ChanSpec{T}(s.corr, filter(kv -> f(kv[1]), s.blocks), s.Δmax, s.chan)
end

Base.length(s::ChanSpec) = length(s.blocks)
Base.size(s::ChanSpec) = size(s.blocks)
nb_blocks(s::ChanSpec) = sum(length(b) for b in s.blocks)

function Base.show(io::IO, s::ChanSpec)
    println(io, "channel $(s.chan), Δmax = $(s.Δmax)")
    nondiags =
        sort([indices(V) for V in fields(s) if !V.diagonal], by = x -> (x[1], x[2]))
    diags = sort(
        [V for V in fields(s) if V.diagonal && !V.degenerate],
        by = V -> real(total_dimension(V)),
    )
    degs = sort(
        [V for V in fields(s) if V.degenerate],
        by = V -> real(total_dimension(V)),
    )
    if !isempty(degs)
        println(io, "Degenerate:")
        for V in degs
            println(io, V)
        end
    end
    if !isempty(diags)
        println(io, "Diagonal:")
        for V in diags
            println(io, V)
        end
    end
    if !isempty(nondiags)
        print(io, "Non-diagonal:")
        r = -Inf
        for rs in nondiags
            if rs[1] != r
                print(io, "\n")
            else
                print(io, ", ")
            end
            r = rs[1]
            print(io, rs)
        end
    end
    println(io, "")
end
