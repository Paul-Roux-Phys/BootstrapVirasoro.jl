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
const Spec = Spectrum

struct BulkSpectrum{T} <: Spectrum{T}
    Δmax::T
    fields::Set{Field{T}}
end
const BulkSpec = BulkSpectrum

struct BdySpectrum{T} <: Spectrum{T}
    Δmax::T
    fields::Set{CD{T}}
end
const BdySpec = BdySpectrum

Spectrum{T}(Δmax::T, ds::Set{CD{T}}) where {T} = BdySpectrum{T}(Δmax, ds)
Spectrum{T}(Δmax::T, ds::Set{Field{T}}) where {T} = BulkSpectrum{T}(Δmax, ds)

function hasdiagonals(s::BulkSpectrum)
    for V in s.fields
        V.diagonal && return true
    end
    return false
end

"""
        Holds all the data for evaluating blocks in a channel.
"""
struct ChannelSpectrum{T}
    Δmax::T
    corr::Corr
    chan::Symbol
    blocks::Dict{Field{T},Block{T}}
end
const ChanSpec = ChannelSpectrum

fields(s::Spectrum) = s.fields
fields(s::ChannelSpectrum) = keys(s.blocks)

function _add_one!(s::BulkSpec, V)
    if !(V in s.fields) && real(total_dimension(V)) < real(s.Δmax)
        push!(s.fields, V)
    end
    return V
end

function _add_one!(s::ChanSpec, V; kwargs...)
    if !(V in fields(s)) && real(total_dimension(V)) < real(s.Δmax)
        b = Block(s.corr, s.chan, V; kwargs...)
        (s.blocks)[V] = b
        return b
    end
end

_add!(s::Union{Spectrum,ChannelSpectrum}, V; kwargs...) = _add_one!(s, V; kwargs...)
_add!(s::Union{Spectrum,ChannelSpectrum}, fields::Vector; kwargs...) =
    foreach(V -> _add_one!(s, V; kwargs...), fields)
_add!(s::Union{Spectrum,ChannelSpectrum}, fields...; kwargs...) =
    foreach(V -> _add_one!(s, V; kwargs...), fields)

"""
        add!(s, fields)
dd one or several fields to the `Spectrum` or `ChannelSpectrum` s, in-place.
"""
add!(s, fields; kwargs...) = _add!(s, fields; kwargs...)

"""
        add(s, fields; kwargs...)
create a new spectrum by copying s and adding one or several fields
to the `Spectrum` or `ChannelSpectrum` s, in-place.
"""
function add(s, fields; kwargs...)
    s2 = deepcopy(s)
    add!(s2, fields; kwargs...)
    return s2
end

function hasdiagonals(s::ChannelSpectrum)
    for V in keys(s.blocks)
        V.diagonal && return true
    end
    return false
end

function Spectrum(
    fields::AbstractArray{Field{T}},
    Δmax;
    interchiral = false,
) where {T}
    s = BulkSpectrum{T}(Δmax, Set{Field{T}}())
    if interchiral
        foreach(V -> if -1 < real(V.s) <= 1 || (!V.isdegenerate && V.diagonal)
            add!(s, V)
        end, fields)
    else
        foreach(V -> add!(s, V), fields)
    end
    return s
end

function Spectrum(fields::AbstractArray{CD{T}}, Δmax; interchiral = false) where {T}
    s = BdySpectrum{T}(Δmax, [])
    if interchiral
        foreach(V -> if -1 < real(V.s) <= 1
            add!(s, V)
        end, fields)
    else
        foreach(V -> add!(s, V), fields)
    end
    return s
end

function ChannelSpectrum(
    co::Co{T},
    Δmax::Number,
    chan;
    kwargs...,
) where {T}
    return ChannelSpectrum{T}(Δmax, co, chan, Dict{Field{T},Block{T}}())
end

function ChannelSpectrum(co::Co, s::Spec{T}, chan; kwargs...) where {T}
    schan = ChannelSpectrum(co, co.Nmax, chan; kwargs...)
    for V in fields(s)
        add!(schan, V; kwargs...)
    end
    return schan
end

function Spectrum(s::ChannelSpectrum{T}) where {T}
    return Spectrum{T}(s.Δmax, Set(fields(s)))
end

function _remove_one!(s::Spectrum, V)
    delete!(s.fields, V)
end

function _remove_one!(s::ChannelSpectrum, V)
    delete!(s.blocks, V)
end

_remove!(s::Union{Spectrum,ChannelSpectrum}, V) = _remove_one!(s, V)
_remove!(s::Union{Spectrum,ChannelSpectrum}, fields::Vector) =
    foreach(V -> _remove_one!(s, V), fields)
_remove!(s::Union{Spectrum,ChannelSpectrum}, fields...) =
    foreach(V -> _remove_one!(s, V), fields)

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

Base.length(s::Spectrum) = length(s.fields)
Base.length(s::ChannelSpectrum) = length(s.blocks)
Base.size(s::Spectrum) = size(s.fields)
Base.size(s::ChannelSpectrum) = size(s.blocks)
nb_blocks(s::ChannelSpectrum) = sum(length(b) for b in s.blocks)

function field_compatible_with_signature(co::Correlation4, V, sig, chan)
    V1, V2, _, _ = permute_4(co.fields, chan)
    return V.diagonal ? (V1.r + V2.r) % 1 == 0 && sig[chan] == 0 :
           (V.r + V1.r + V2.r) % 1 == 0 && V.r >= sig[chan]
end

function field_compatible_with_signature(co::Correlation1, V, sig, chan)
    return V.diagonal ? sig[chan] == 0 :
           V.r % 1 == sig[chan] % 1 && V.r >= sig[chan]
end

function ChannelSpectra(
    co,
    s::Channels{<:Spectrum{T}},
    signature = (s = 0, t = 0, u = 0);
    kwargs...,
) where {T}
    chans = (:s, :t, :u)
    schan = (;
        (chan => ChannelSpectrum(co, s.s.Δmax, chan; kwargs...) for chan in chans)...
    )
    exclude = nothing
    if haskey(kwargs, :exclude)
        exclude = kwargs[:exclude]
    end
    kwargs = (; (k => v for (k, v) in kwargs if k != :exclude)...)

    # Launch one task per channel
    tasks = map(chans) do chan
        Threads.@spawn begin
            for V in s[chan].fields
                cond = field_compatible_with_signature(co, V, signature, chan)
                if haskey(kwargs, :parity) && kwargs[:parity] != 0
                    cond = cond && (real(V.s) >= 0 || V.diagonal)
                end
                if !isnothing(exclude) && haskey(exclude, chan)
                    cond = cond && !(V in exclude[chan])
                end
                if cond
                    add!(schan[chan], V; kwargs...)
                end
            end
        end
    end

    # Wait for all threads to finish
    foreach(fetch, tasks)

    return schan
end

ChannelSpectra(co, S::Spectrum, signature = (s = 0, t = 0, u = 0); kwargs...) =
    ChannelSpectra(co, (s = S, t = S, u = S), signature; kwargs...)

function Base.show(io::IO, s::BulkSpectrum)
    nondiags =
        sort([indices(V) for V in s.fields if !V.diagonal], by = x -> (x[1], x[2]))
    diags = sort(
        [V for V in s.fields if V.diagonal && !V.isdegenerate],
        by = V -> real(total_dimension(V)),
    )
    degs = sort(
        [V for V in s.fields if V.isdegenerate],
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

function Base.show(io::IO, s::ChanSpec)
    println(io, "channel $(s.chan), Δmax = $(s.Δmax)")
    show(io, Spectrum(s))
end
