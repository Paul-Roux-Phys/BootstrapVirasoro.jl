function Cref(V₁, V₂, V₃, DG)
    β = V₁.c.β
    r₁, s₁ = get_indices(V₁)
    r₂, s₂ = get_indices(V₂)
    r₃, s₃ = get_indices(V₃)

    return prod(
        1 / DG(
            (β + 1 / β) / 2 +
            β / 2 * abs(pm₁ * r₁ + pm₂ * r₂ + pm₃ * r₃) +
            1 / 2 / β * (pm₁ * s₁ + pm₂ * s₂ + pm₃ * s₃),
        ) for pm₁ in (-1, 1) for pm₂ in (-1, 1) for pm₃ in (-1, 1)
    )
end

Cref(V₁, V₂, V₃) = Cref(V₁, V₂, V₃, DoubleGamma(V₁.c.β))

function Bref(DG, c, r, s, reg = 1 / big(10^15))
    β = c.β
    if r % 1 == 0 && s % 1 == 0
        s += reg
    end
    π = oftype(c.β, Base.π) # π in the correct precision
    return (-1)^(round(Int, r * s)) / 2 / sin(π * (r % 1 + s)) /
           sin(π * (r + s / β^2)) / prod(
        DG(β + pm1 * β * r + pm2 * s / β) for pm1 in (-1, 1) for pm2 in (-1, 1)
    )
end

function Bref(DG, c, P)
    β = c.β
    prod(1 / DG(β^pm1 + pm2 * 2P) for pm1 in (-1, 1) for pm2 in (-1, 1))
end

function Bref(V::Field, DG, reg = 1 / big"10"^15)
    c = V.c
    if V.diagonal
        return Bref(DG, c, V.dims[:left].P)
    else
        return Bref(DG, c, V.r, V.s, reg)
    end
end

function compute_reference(co::Correlation4, V::Field, chan, DG)
    V₁, V₂, V₃, V₄ = getfields(co, chan)
    return Cref(V₁, V₂, V, DG) * Cref(V₃, V₄, V, DG) / Bref(V, DG)
end

function compute_reference(co::Correlation1, V::Field, chan, DG)
    V₁ = getfields(co, chan)
    return Cref(V₁, V, V, DG) / Bref(V, DG)
end

compute_reference(b::Block, DG) = compute_reference(b.corr, b.chan_field, b.chan, DG)

struct StructureConstants{T<:Complex}
    constants::Channels{Dict{Field{T},T}}
    errors::Channels{Dict{Field{T},Float32}}
    reference::Channels{Dict{Field{T},T}}
end
const StrCst = StructureConstants
const SC = StructureConstants

function StructureConstants{T}() where {T}
    constants = @channels Dict{Field{T}, T}()
    errors = @channels Dict{Field{T}, Float32}()
    reference = deepcopy(constants)
    return StructureConstants{T}(constants, errors, reference)
end

function Base.getproperty(c::SC, s::Symbol)
    s === :fields && begin
        consts = getfield(c, :constants)
        return vcat([[V for V in keys(consts[chan])] for chan in keys(consts)]...)
    end
    getfield(c, s)
end

Base.getindex(c::SC, s::Symbol) = c.constants[s]

function fix!(cnst, chan, field, value; error = 0)
    cnst[chan][field] = value
    cnst.errors[chan][field] = error
end

Base.length(c::SC) =
    sum(length(c.constants[chan]) for chan in keys(c.constants))

function compute_reference!(c::SC, b::Block, DG)
    chan = b.channel
    V = b.channel_field
    c.reference[chan][V] = compute_reference(b, DG)
end

function compute_reference!(c::SC, S::ChanSpec, DG)
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

function compute_reference!(c::SC, S::Channels{<:ChanSpec})
    β = S[1].corr.c.β
    DG = DoubleGamma(β)
    for chan in keys(S)
        compute_reference!(c, S[chan], DG)
    end
    return c.reference
end

function to_dataframe(c::SC, rmax = nothing)
    df = DataFrame(
        Channel = Symbol[],
        Field = String[],
        StructureConstant = Complex[],
        RelativeError = Float32[],
    )

    for chan in CHANNELS
        fields =
            sort(collect(keys(c.constants[chan])), by = V -> real(total_dimension(V)))
        rmax !== nothing && (fields = filter(V -> V.r <= rmax, fields))
        for V in fields
            push!(df, (chan, string(V), c.constants[chan][V], c.errors[chan][V]))
        end
    end

    return df
end

function Base.show(io::IO, c::SC; rmax = nothing, plaintext = false)
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
        channel_header = " Channel $(table.Channel[1])\n"
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
            printstyled(io, channel_header; bold = true)
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

function save(io::IO, c::SC; format::Symbol = :csv)
    df = to_dataframe(c)
    if format == :csv
        CSV.write(io, df)
    else
        show(io, c, plaintext = true)  # fallback to default text show
    end
end
