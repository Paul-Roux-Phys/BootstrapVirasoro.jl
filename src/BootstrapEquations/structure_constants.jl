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
    constants = (; (chan => Dict() for chan in (:s, :t, :u))...)
    errors = (; (chan => Dict() for chan in (:s, :t, :u))...)
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
           real(total_dimension(V)) < real(S.Δmax) / 3
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

function to_dataframe(c::StructureConstants, rmax = nothing)
    df = DataFrame(
        Channel = Symbol[],
        Field = String[],
        StructureConstant = Complex[],
        RelativeError = Float32[],
    )

    for chan in sort(collect(keys(c.constants)), by = string)
        fields =
            sort(collect(keys(c.constants[chan])), by = V -> real(total_dimension(V)))
        rmax !== nothing && (fields = filter(V -> V.r <= rmax, fields))
        for V in fields
            push!(df, (chan, string(V), c.constants[chan][V], c.errors[chan][V]))
        end
    end

    return df
end

function Base.show(io::IO, c::StructureConstants; rmax = nothing, plaintext = false)
    df = to_dataframe(c, rmax)
    grouped = groupby(df, :Channel)
    for table in grouped
        t = Tables.columntable(select(table, Not(:Channel)))
        hl_conv = Highlighter(
            (table, i, j) ->
                table[:RelativeError][i] < 2^(-precision(BigFloat, base = 10) / 2) &&
                    j >= 2,
            crayon"green",
        )
        hl_zero = Highlighter(
            (table, i, j) ->
                table[:RelativeError][i] > 1e-2 &&
                    abs(table[:StructureConstant][i]) < 2^(-precision(BigFloat, base = 10) / 3) &&
                    j >= 2,
            crayon"green",
        )
        hl_not_conv = Highlighter(
            (table, i, j) ->
                table[:RelativeError][i] > 1e-2 &&
                    abs(table[:StructureConstant][i]) > 2^(-precision(BigFloat, base = 10) / 3) &&
                    j >= 2,
            crayon"red",
        )
        fmt_zero =
            (v, i, j) ->
                (
                    t.RelativeError[i] > 1e-2 &&
                    abs(t.StructureConstant[i]) < 2^(-precision(BigFloat, base = 10) / 3) &&
                    j == 2
                ) ? 0 : v
        channel_header = "="^20 * "\nChannel $(table.Channel[1])\n" * "="^20 * "\n"
        if plaintext
            print(io, channel_header)
            pretty_table(
                io,
                t,
                header = ["Field", "Structure constant", "Relative error"],
                header_alignment = :l,
                formatters = (fmt_zero,),
            )
        else
            printstyled(io, channel_header, color = :red, bold = true)
            pretty_table(
                io,
                t,
                header = ["Field", "Structure constant", "Relative error"],
                header_alignment = :l,
                header_crayon = crayon"blue",
                highlighters = (hl_conv, hl_zero, hl_not_conv),
                formatters = (fmt_zero,),
            )
        end
    end
end

function save(io::IO, c::StructureConstants; format::Symbol = :csv)
    df = to_dataframe(c)
    if format == :csv
        CSV.write(io, df)
    else
        show(io, c, plaintext = true)  # fallback to default text show
    end
end
