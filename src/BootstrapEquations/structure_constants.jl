include("reference_constants.jl")

"""
TODO
"""
struct StructureConstants{T<:Complex}
    constants::Channels{Dict{Field{T},T}}
    errors::Channels{Dict{Field{T},Float32}}
    reference::Channels{Dict{Field{T},T}}
end

function StructureConstants{T}() where {T}
    constants = Channels{Dict{Field{T},T}}(Tuple(Dict() for chan in (:s, :t, :u)))
    errors = Channels{Dict{Field{T},Float32}}(Tuple(Dict() for chan in (:s, :t, :u)))
    reference = deepcopy(constants)
    return StructureConstants{T}(constants, errors, reference)
end

function Base.getproperty(c::StructureConstants, s::Symbol)
    s === :fields && begin
        consts = getfield(c, :constants)
        return vcat([[V for V in keys(consts[chan])] for chan in keys(consts)]...)
    end
    getfield(c, s)
end

Base.getindex(c::StructureConstants, s::Symbol) = c.constants[s]

function fix!(cnst, chan, field, value; error = 0)
    cnst[chan][field] = value
    cnst.errors[chan][field] = error
end

Base.length(c::StructureConstants) =
    sum(length(c.constants[chan]) for chan in keys(c.constants))

function compute_reference!(c::StructureConstants, b::Block, DG)
    chan = b.channel
    V = b.channel_field
    c.reference[chan][V] = compute_reference(b, DG)
end

function compute_reference!(c::StructureConstants, S::ChannelSpectrum, DG)
    for b in values(S.blocks)
        chan = b.channel
        V = b.channel_field
        if V in keys(c.constants[chan]) &&
           c.errors[chan][V] < 1e-4 &&
           real(total_dimension(V)) < real(S.Δmax)/3
            compute_reference!(c, b, DG)
        end
    end
end

function compute_reference!(c::StructureConstants, S::Channels{<:ChannelSpectrum})
    β = S[1].corr.c.β
    DG = DoubleGamma(β)
    for chan in keys(S)
        compute_reference!(c, S[chan], DG)
    end
    return c.reference
end

function to_dataframe(c::StructureConstants)
    df = DataFrame(
        Channel = Symbol[],
        Field = String[],
        StructureConstant = Complex[],
        RelativeError = Float32[],
    )

    for chan in sort(collect(keys(c.constants)), by = string)
        fields = collect(keys(c.constants[chan]))
        for V in sort(fields, by = V->real(total_dimension(V)))
            push!(df, (chan, string(V), c.constants[chan][V], c.errors[chan][V]))
        end
    end

    return df
end

function Base.show(io::IO, c::StructureConstants; chan = nothing)
    df = to_dataframe(c)
    show(io, df)
end

function save(io::IO, c::StructureConstants; format::Symbol = :csv)
    df = to_dataframe(c)
    if format == :csv
        CSV.write(io, df)
    else
        show(io, c)  # fallback to default text show
    end
end
