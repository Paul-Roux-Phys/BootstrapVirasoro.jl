export Spectrum,
    ChannelSpectrum,
    add!, remove!,
    nb_blocks

"""
    Spectrum{T}

Abstract type for representing a CFT spectrum.
"""
abstract type Spectrum{T} end

mutable struct BulkSpectrum{T} <: Spectrum{T}
    Δmax::T
    interchiral::Bool
    fields::Vector{Field{T}}
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
end

function _add_one!(s::ChannelSpectrum, V)
    if !(V in s.fields) && real(total_dimension(V)) < real(s.Δmax)
        push!(s.blocks, Block(s.corr, s.channel, V; interchiral=s.interchiral, Δmax=s.Δmax))
    end
end

add!(s::Union{Spectrum, ChannelSpectrum}, V) = _add_one!(s, V)
add!(s::Union{Spectrum, ChannelSpectrum}, fields::Vector) = foreach(V -> _add_one!(s, V), fields)
add!(s::Union{Spectrum, ChannelSpectrum}, fields...) = foreach(V -> _add_one!(s, V), fields)

function Spectrum(fields::Vector{Field{T}}, Δmax; interchiral=false) where {T}
    s = BulkSpectrum{T}(Δmax, interchiral, [])
    foreach(V -> add!(s, V), fields)
    return s
end

function ChannelSpectrum(co::Correlation, s::Spectrum{T}, chan) where {T}
    schan = ChannelSpectrum{T}(s.Δmax, co, chan, s.interchiral, [])
    foreach(V -> add!(schan, V), s.fields)
    return schan
end

function ChannelSpectra(co, s::Spectrum{T}; signature=Dict(:s => 0, :t => 0, :u => 0)) where {T}
    schan = Dict(
        chan => ChannelSpectrum{T}(s.Δmax, co, chan, s.interchiral, [])
        for chan in (:s, :t, :u)
    )
    for V in s.fields
        for chan in (:s, :t, :u)
            if isdegenerate(V) || V.r >= signature[chan]
                add!(schan[chan], V)
            end
        end
    end
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
end

function _remove_one!(s::ChannelSpectrum, V)
    idx = findfirst(==(V), s.fields)
    if idx !== nothing
        deleteat!(s.blocks, idx)
    end
end

remove!(s::Union{Spectrum, ChannelSpectrum}, V) = _remove_one!(s, V)
remove!(s::Union{Spectrum, ChannelSpectrum}, fields::Vector) = foreach(V -> _remove_one!(s, V), fields)
remove!(s::Union{Spectrum, ChannelSpectrum}, fields...) = foreach(V -> _remove_one!(s, V), fields)

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
