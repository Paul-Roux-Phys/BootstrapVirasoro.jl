export Spectra,
    Spectrum,
    BulkSpectrum,
    blocks, nb_blocks,
    generate_pairs

"""
    Spectrum{T}

Abstract type for representing a CFT spectrum.
"""
abstract type Spectrum{T} end

const Spectra{T} = Dict{Symbol, Spectrum{T}}

mutable struct BulkSpectrum{T} <: Spectrum{T}
    
    fields::Vector{Field{T}}
    blocks::Vector{Block{T}}

end

function blocks(co::Correlation, chan, Vs::Vector{Field{T}}; kwargs...) where {T}
    [Block(co, chan, V; kwargs...) for V in Vs]
end

function Spectrum(co, chan, Vs::Vector{Field{T}}; kwargs...) where {T}
    bs = blocks(co, chan, Vs; kwargs...)
    BulkSpectrum{T}(Vs, bs)
end

function Base.getproperty(s::Spectrum, p::Symbol)
    p === :chan && return getfield(s, blocks)[1].channel
    getfield(s, p)
end

function add!(s::BulkSpectrum, V::Field)
    if !(V in s.fields)
        push!(s.fields, V)
    end
end

function add!(s::BulkSpectrum, fields::Vector{Field})
    for V in fields
        add!(s, V)
    end
end

function add!(s::BulkSpectrum, Vs::Field...)
    for V in Vs
        add!(s, V)
    end
end

Base.length(s::Spectrum) = length(s.fields)
Base.size(s::Spectrum) = size(s.fields)
nb_blocks(s::Spectrum) = sum(length(b) for b in s.blocks)

function Base.show(io::IO, s::BulkSpectrum)
    nondiags = sort(
        [V.indices for V in s.fields if V.r != 0],
        by = x -> (x[1], x[2])
    )
    diags = sort(
        [V for V in s.fields if V.r == 0],
        by = abs
    )
    if !isempty(diags)
        println(io, "Diagonal part:")
        for V in diags
            println(io, V)
        end
    end
    if !isempty(nondiags)
        print(io, "Non-diagonal part indices:")
        r = -Inf
        for rs in nondiags
            if rs[1] != r
                println(io, "")
            else
                print(io, ", ")
            end
            r = rs[1]
            print(io, rs)
        end
    end
end

function generate_pairs(r_range, conditions=(r, s)->(r*s%1==0))
    [(r, s) for r in r_range for s in -1+1//r:1//r:1 if conditions(r, s)]
end

# condition = (r, s) -> (r%1 == 0 && r*s%1 == 0)

# generate_pairs(1//2:1//2:3, -3:1//2:3, (r, s) -> ((r*s)%1 == 0 && r%1 == 0))
