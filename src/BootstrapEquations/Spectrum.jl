struct ChannelSpectrum{T}
    Δmax::T
    corr::Corr
    chan::Symbol
    blocks::Dict{Field{T},Block{T}}
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

function ChannelSpectrum(
    co::Co{T},
    chan,
    Vs::AbstractArray{<:Field},
    f::Function,
) where {T}
    s = ChannelSpectrum{T}(co.Δmax, co, chan, Dict{Field{T},Block{T}}())
    Threads.@threads for V in Vs
        add!(s, f(V))
    end
    return s
end

function _remove_one!(s::ChanSpec, V)
    delete!(s.blocks, V)
end

_remove!(s::ChanSpec, V) = _remove_one!(s, V)
_remove!(s::ChanSpec, fields::Vector) = foreach(V -> _remove_one!(s, V), fields)
remove!(s, V) = _remove!(s, V)
remove(s, Vs) = remove!(deepcopy(s), Vs)

Base.length(s::ChanSpec) = length(s.blocks)
Base.size(s::ChanSpec) = size(s.blocks)
nb_blocks(s::ChanSpec) = sum(length(b) for b in s.blocks)

function LoopBlock(co, chan, V, interchiral, parity)
    if parity == 0 || V.diagonal
        return Block(co, chan, V, interchiral = interchiral)
    elseif parity == 1 || parity == -1 && V.s > 0
        return Block(co, chan, V, interchiral = interchiral) +
               parity * Block(co, chan, reflect(V), interchiral = interchiral)
    end
end

function parity_compat(co::Correlation4, V::Field, chan)
    V1, V2, _, _ = permute_4(co.fields, chan)
    r = get_indices(V)[1]
    return (r + V1.r + V2.r) % 1 == 0
end

parity_compat(_::Correlation1, _::Field, _) = true

function LoopSpectrum(co, chan, Vs, sig, interchiral = true, parity = 0)
    Vs = filter(
        V -> real(total_dimension(V)) < co.Δmax && get_indices(V)[1] >= sig &&
            parity_compat(co, V, chan),
        Vs
    )
    if parity != 0
        Vs = filter(V -> V.diagonal || V.s >= 0, Vs)
    end
    LB = (V, chan) -> LoopBlock(co, chan, V, interchiral, parity)
    return ChanSpec(co, chan, Vs, V -> LB(V, chan))
end

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
