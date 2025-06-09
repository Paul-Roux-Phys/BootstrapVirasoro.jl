"""
    Spectrum(fields, Δmax; interchiral=false)

Abstract type for representing a CFT spectrum.

# examples

```jldoctest
julia> c = CentralCharge(β=1/(big"0.8"+big"0.1"*im));

julia> V1 = Field(c, r=2, s=0);

julia> co = Correlation(V1, Δmax=10.);

julia> fields = [Field(c, r=r, s=s) for r in 2:30 for s in -1+1//r:1//r:1 if r*s% 1 == 0];

julia> s = Spectrum(fields, 10.0)
Non-diagonal:
(2, -1//2), (2, 0), (2, 1//2), (2, 1)
(3, -2//3), (3, -1//3), (3, 0), (3, 1//3), (3, 2//3), (3, 1)

julia> add!(s, Field(c, r=1, s=0)); s
Non-diagonal:
(1, 0)
(2, -1//2), (2, 0), (2, 1//2), (2, 1)
(3, -2//3), (3, -1//3), (3, 0), (3, 1//3), (3, 2//3), (3, 1)
```
"""
abstract type Spectrum{T} end

mutable struct BulkSpectrum{T} <: Spectrum{T}
    Δmax::T
    interchiral::Bool
    fields::Vector{Field{T}}
end

mutable struct BoundarySpectrum{T} <: Spectrum{T}
    Δmax::T
    interchiral::Bool
    fields::Vector{ConformalDimension{T}}
end

"""
        Holds all the data for evaluating blocks in a channel.
"""
mutable struct ChannelSpectrum{T}
    Δmax::T
    corr::Correlation{T}
    channel::Symbol
    interchiral::Bool
    blocks::Vector{Block{T}}
end

function _add_one!(s::Spectrum, V)
    if !(V in s.fields) && real(total_dimension(V)) < real(s.Δmax)
        push!(s.fields, V)
    end
    return V
end

function _add_one!(s::ChannelSpectrum, V)
    if !(V in s.fields) && real(total_dimension(V)) < real(s.Δmax)
        b = Block(s.corr, s.channel, V; interchiral=s.interchiral, Δmax=s.Δmax)
        push!(s.blocks, b)
        return b
    end
end

_add!(s::Union{Spectrum, ChannelSpectrum}, V) = _add_one!(s, V)
_add!(s::Union{Spectrum, ChannelSpectrum}, fields::Vector) = foreach(V -> _add_one!(s, V), fields)
_add!(s::Union{Spectrum, ChannelSpectrum}, fields...) = foreach(V -> _add_one!(s, V), fields)

"""
        add!(s, fields)
add one or several fields to the `Spectrum` or `ChannelSpectrum` s, in-place.
"""
add!(s, fields) = _add!(s, fields)

"""
        add(s, fields)
create a new spectrum by copying s and adding one or several fields
to the `Spectrum` or `ChannelSpectrum` s, in-place.
"""
function add(s, fields)
    s2 = deepcopy(s)
    add!(s2, fields)
    return s2
end

function Spectrum(fields::Vector{Field{T}}, Δmax; interchiral=false) where {T}
    s = BulkSpectrum{T}(Δmax, interchiral, [])
    if interchiral
        foreach(V -> if -1 < V.s <= 1 add!(s, V) end, fields)
    else
        foreach(V -> add!(s, V), fields)
    end
    return s
end

function Spectrum(fields::Vector{ConformalDimension{T}}, Δmax; interchiral=false) where {T}
    s = BoundarySpectrum{T}(Δmax, interchiral, [])
    foreach(V -> add!(s, V), fields)
    return s
end

function ChannelSpectrum(co::Correlation, s::Spectrum{T}, chan) where {T}
    schan = ChannelSpectrum{T}(s.Δmax, co, chan, s.interchiral, [])
    for i in eachindex(s.fields)
        add!(schan, s.fields[i])
    end
    return schan
end

function ChannelSpectra(co, s::Spectrum{T}; signature=nothing) where {T}
    chans = (:s, :t, :u)
    schan = Channels{ChannelSpectrum{T}}(Tuple(
        ChannelSpectrum{T}(s.Δmax, co, chan, s.interchiral, [])
        for chan in chans
    ))

    # Launch one task per channel
    tasks = map(chans) do chan
        Threads.@spawn begin
            V1, V2, V3, V4 = permute_fields(co.fields, chan)
            for V in s.fields
                cond = isdiagonal(V) ? (V1.r + V2.r) % 1 == 0 :
                                       (V.r + V1.r + V2.r) % 1 == 0
                if cond && V.r >= signature[chan]
                    add!(schan[chan], V)
                end
            end
        end
    end

    # Wait for all threads to finish
    foreach(fetch, tasks)

    return schan
end

function Spectrum(s::ChannelSpectrum{T}) where {T}
    return Spectrum(Vector{Field{T}}(s.fields), s.Δmax)
end

function Base.getproperty(s::ChannelSpectrum, p::Symbol)
    p === :fields && return [b.channel_field for b in s.blocks]
    getfield(s, p)
end

function _remove_one!(s::Spectrum, V)
    idx = findfirst(==(V), s.fields)
    if idx !== nothing
        deleteat!(s.fields, idx)
    end
    return idx
end

function _remove_one!(s::ChannelSpectrum, V)
    idx = findfirst(==(V), s.fields)
    if idx !== nothing
        deleteat!(s.blocks, idx)
    end
    return idx
end

_remove!(s::Union{Spectrum, ChannelSpectrum}, V) = _remove_one!(s, V)
_remove!(s::Union{Spectrum, ChannelSpectrum}, fields::Vector) = foreach(V -> _remove_one!(s, V), fields)
_remove!(s::Union{Spectrum, ChannelSpectrum}, fields...) = foreach(V -> _remove_one!(s, V), fields)

"""
        remove!(s, V)
remove one or several fields from the `Spectrum` or `ChannelSpectrum` s.
"""
remove!(s, V) = _remove!(s, V)

"""
        remove(s, fields)
create a new spectrum by copying and removing one or several fields to the `Spectrum` or `ChannelSpectrum` s.
"""
function remove(s, fields)
    s2 = deepcopy(s)
    remove!(s2, fields)
    return s2
end


Base.length(s::Union{Spectrum, ChannelSpectrum}) = length(s.fields)
Base.size(s::Union{Spectrum, ChannelSpectrum}) = size(s.fields)
nb_blocks(s::ChannelSpectrum) = sum(length(b) for b in s.blocks)

function Base.show(io::IO, s::BulkSpectrum)
    nondiags = sort(
        [V.indices for V in s.fields if !isdiagonal(V)],
        by = x -> (x[1], x[2])
    )
    diags = sort(
        [V for V in s.fields if isdiagonal(V) && !isdegenerate(V)],
        by = V -> real(total_dimension(V))
    )
    degs =  sort(
        [V for V in s.fields if isdegenerate(V)],
        by = V -> real(total_dimension(V))
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
end

function Base.show(io::IO, s::ChannelSpectrum)
    println(io, "channel $(s.channel), Δmax = $(s.Δmax)")
    show(io, Spectrum(s))
end
