export Spectrum,
    BulkSpectrum

"""
    Spectrum{T}

Abstract type for representing a CFT spectrum.
"""
abstract type Spectrum{T} end

mutable struct BulkSpectrum{T} <: Spectrum{T}
    
    channel::Symbol
    maxdim::ConformalDimension{T}
    fields::Vector{Field{T}}
    blocks::Dict{Field{T}, GeneralizedBlock{T}}

end

function add!(s::BulkSpectrum, V::Field; interchiral=false)
    if !(V in s.fields)
        push!(s.fields, V)
        s.blocks[V] = GeneralizedBlock(c, channel, field, Δmax, interchiral=interchiral)
    end
end

function add!(s::BulkSpectrum, fields::Vector{Field}; interchiral=false)
    for V in fields
        add!(s, V, interchiral=interchiral)
    end
end

function add!(s::BulkSpectrum, Vs::Field...; interchiral=false)
    for V in Vs
        add!(s, V, interchiral=interchiral)
    end
end

function BulkSpectrum(
    c::Correlation{T},
    channel::Symbol,
    Δmax::ConformalDimension,
    fields;
    interchiral=false
) where {T}
    res = BulkSpectrum{T}(channel, Δmax, [], Dict())

    for field in fields
        res.blocks[field] = GeneralizedBlock(c, channel, field, Δmax, interchiral=interchiral)
    end
    res.fields = fields

    res
end

function Base.show(io::IO, s::BulkSpectrum)
    nondiags = sort(
        [V.indices for V in s.fields if !V.isdiagonal],
        by = x -> (x[1], x[2])
    )
    diags = sort(
        [V for V in s.fields if V.isdiagonal],
        by = abs
    )
    println(io, "Diagonal part:")
    for V in diags
        println(io, V)
    end
    println(io, "Non-diagonal-part:")
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

function generate_pairs(r_range, s_range, condition)
    [(r, s) for r in r_range for s in -1:1//r:1 if condition(r, s)]
end

condition = (r, s) -> (r%1 == 0 && r*s%1 == 0)

generate_pairs(1//2:1//2:3, -3:1//2:3, (r, s) -> ((r*s)%1 == 0 && r%1 == 0))