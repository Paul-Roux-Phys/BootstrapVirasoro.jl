"""
TODO
"""
struct StructureConstants{T}
    constants::Channels{Dict{Field{T},T}}
    errors::Channels{Dict{Field{T},T}}
end

function StructureConstants{T}() where {T}
    constants = Channels{Dict{Field{T},T}}(Tuple(Dict() for chan in (:s, :t, :u)))
    errors = deepcopy(constants)
    return StructureConstants{T}(constants, errors)
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

function Base.show(io::IO, c::StructureConstants{T}) where {T}
    chans = keys(c.constants)
    nondiags = Dict(
        chan =>
            sort([V for V in keys(c[chan]) if !isdiagonal(V)], by = V -> (V.r, V.s)) for
        chan in chans
    )
    diags = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if isdiagonal(V) && !isdegenerate(V)],
            by = V -> real(total_dimension(V)),
        ) for chan in chans
    )
    degs = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if isdegenerate(V)],
            by = V -> real(total_dimension(V)),
        ) for chan in chans
    )

    for chan in sort(collect(chans), by = string)
        # Collect all the labels for this channel
        all_Vs = vcat(degs[chan], diags[chan], nondiags[chan])

        # Find max width of label for alignment
        max_label_width = if isempty(all_Vs)
            0
        else
            maximum(length(string(V)) for V in all_Vs)
        end

        str_cst_col_width = if isempty(all_Vs)
            0
        else
            length(format_complex(zero(T)))
        end

        # Print channel header
        println(io, repeat('=', 9 + length(string(chan))))
        println(io, "Channel $(chan)")
        println(io, repeat('=', 9 + length(string(chan))))

        # Print column headers
        label1 = rpad("Field", max_label_width)
        label2 = rpad("Structure constant", str_cst_col_width)
        label3 = rpad("Relative error", str_cst_col_width)
        println(io, "$label1 | $label2  | $label3")
        println(io, repeat("-", max_label_width+2*str_cst_col_width))

        # Print each row with alignment
        for V in all_Vs
            label = rpad(string(V), max_label_width)
            value = format_complex(c[chan][V])
            error = @sprintf("%.2e", c.errors[chan][V])
            println(io, "$label | $value | $error")
        end
    end
end

function write_csv(io::IO, c::StructureConstants{T}) where {T}
    chans = keys(c.constants)
    cc = collect(keys(c[:s]))[1].c
    println("centralcharge=$(cc)")
    println(
        io,
        "chan,r,s,isdiagonal,isdegenerate,real(structureconstant)," *
        "imag(structureconstant)," *
        "err",
    )
    for chan in sort(collect(chans), by = string)
        for V in keys(c[chan])
            if !haskey(c.errors[chan], V)
                continue
            end
            sc = c[chan][V]
            err = real(c.errors[chan][V])
            println(
                io,
                "$(chan),$(V.r),$(V.s),$(isdiagonal(V)),$(isdegenerate(V))," *
                "$(real(sc)),$(imag(sc))," *
                "$err",
            )
        end
    end
end

Base.length(c::StructureConstants) =
    sum(length(c.constants[chan]) for chan in keys(c.constants))
