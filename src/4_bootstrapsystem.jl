struct FieldCache
    block_values::Vector{Acb}
    factors     ::Vector{Acb}
end

struct StructureConstants
    dict::Dict{CDorField,Tuple{Acb,Float64}}
end
const StructConst = StructureConstants
const SC = StructureConstants

"""
data required to build a system of bootstrap equations

# Fields

unknowns: list of fields whose structure constant is unknown, in a Channels struct
knowns: list of known (fields, structure_constant), in a Channels struct
block_values: values of the conformal blocks, in a Channels struct
factors: prefactors multiplying the block values. May depend on the moduli, the fields
"""
struct BootstrapSystem
    unknowns     ::Channels{Vector{CDorField}}
    knowns       ::Channels{Vector{Tuple{CDorField,Acb}}}
    moduli_cache ::Channels
    fields_cache ::Channels{Dict{CDorField,FieldCache}}
end

function BootstrapSystem(co::Correlation, blocks;
                         fix=nothing, moduli=nothing, extrapoints::Int=6,
                         prefactors=nothing)
    chans = collect(keys(blocks))
    chan_count = length(chans)

    unknowns_set = Channels(Dict(chan => Set{CDorField}(b.chan_field for b in blocks[chan])
                                 for chan in chans))
    knowns = Channels(Dict(chan => Tuple{CDorField, Acb}[] for chan in chans))

    # if nothing is fixed, fix the structure constant of the first found diagonal field to 1
    for chan in chans
        if fix !== nothing
            break
        end
        for V in unknowns_set[chan]
            if V.diagonal
                fix = [(chan, V, one(Acb))]
                break
            end
        end
    end
    # if no diagonal field is found, fix the first structure constant we find in the s-channel to 1
    if fix === nothing
        V = first(unknowns_set.s)
        fix = [(:s, V, one(Acb))]
    end

    for (chan, V, D) in fix
        delete!(unknowns_set[chan], V)
        push!(knowns[chan], (V, D))
    end
    unknowns = Channels(Dict(chan => collect(unknowns_set[chan])
                             for chan in chans))
    
    unknowns_count = sum(length(unknowns[chan]) for chan in chans)
    lines_count = 2 * ((unknowns_count + extrapoints) ÷ 2)
    if moduli !== nothing
        moduli_count = length(moduli)
        if !(moduli_count > lines_count ÷ (chan_count - 1))
            error("You must give enough positions for the system to be overdetermined")
        end
    else
        moduli_count = chan_count == 1 ? lines_count ÷ 2 : lines_count ÷ (chan_count - 1)
    end

    moduli = choose_random_moduli(moduli_count, co)
     
    # compute the cache for all moduli, in each channel
    moduli_chan = Channels(Dict(chan => [channel_modulus(x, co, chan) for x in moduli]
                             for chan in chans))
    moduli_cache = Channels(Dict(chan => Vector(undef, length(moduli_chan[chan]))
                              for chan in chans)) # result container
    Threads.@threads for i = 1:chan_count
        chan = chans[i]
        moduli_cache[chan] = [PosCache(x, co, chan) for x in moduli_chan[chan]]
    end

    # compute the cached data for each field
    fields_cache = Channels(Dict(chan => Dict{CDorField,FieldCache}()
                                 for chan in chans))

    # evaluate the blocks and prefactors
    for chan in chans
        thread_results = [Dict{CDorField,Vector{Acb}}() for _ = 1:max_thread_id()] # thread-local results
        Threads.@threads for i = 1:length(blocks[chan])
            b = blocks[chan][i]
            thread_results[Threads.threadid()][b.chan_field] =
                [b(x) for x in moduli_cache[chan]]
        end
        # Merge all thread-local dictionaries, and compute prefactors
        for d in thread_results
            for (V, vals) in d
                if prefactors === nothing
                    fact = ones(Acb, moduli_count)
                else
                    fact = [prefactors(chan, V, moduli[k]) for k in 1:moduli_count]
                end
                fields_cache[chan][V] = FieldCache(vals, fact)
            end
        end
    end

    return BootstrapSystem(unknowns, knowns, moduli_cache, fields_cache)
end

function new_random_point(; square=true)
    xmin, xmax = square ? (0.1, 0.5) : (-0.4, 1.4)
    return z = xmin + (xmax - xmin) * rand() + (1 + 4 * rand()) * im / 10
end

function new_random_point!(random_points, N; square = true, cond = p -> true)
    while true
        z = new_random_point(square=square)
        sep = square ? 0.2 : 0.4
        if minimum(append!(abs.(z .- random_points), 1)) > sep / sqrt(N) &&
           cond(z) == true
            push!(random_points, z)
            break
        end
    end
end

"""
Generates N points in the square (.1, .5) + (.1, .5)i, while respecting 
a minimum distance .2/sqrt(N) between points. Or in the rectangle 
(-.4, 1.4) + (.1, .5)i with a distance .4/sqrt(N).
N = the number of points.
function = a function that we may apply on the points
square = whether to use a square (otherwise, a rectangle)
"""
function random_points(N; transfo = nothing, square = true, cond = p -> true)
    # Set a temporary seed
    rng = copy(Random.default_rng())
    Random.seed!(3456)
    res = []
    for _ = 1:N
        new_random_point!(res, N, square = square, cond = cond)
    end
    return transfo !== nothing ? transfo.(res) : res

    # Restore RNG state
    Random._GLOBAL_RNG[] = rng
end

choose_random_moduli(moduli_count, _::Correlation1) = 
    Vector{Acb}(random_points(moduli_count, transfo = x -> τfromx(x)))
choose_random_moduli(moduli_count, _::Correlation4) = 
    Vector{Acb}(random_points(moduli_count, transfo = nothing))

function crossratio(chan, x)
    chan === :s && return x
    chan === :t && return 1 - x
    chan === :u && return 1 / x
    error("""Incorrect channel specification in crossratio(channel, x):
          must be in (:s, :t, :u)""")
end

function modular_param(chan, τ)
    chan === :s && return τ
    chan === :t && return -1 / τ
    chan === :u && return (τ - 2) / (τ - 1)
    error("""Incorrect channel specification in crossratio(channel, x):
          must be in (:s, :t, :u)""")
end

channel_modulus(x, _::Correlation4, chan) = crossratio(chan, x)
channel_modulus(x, _::Correlation1, chan) = modular_param(chan, x)

struct BootstrapMatrix
    fields::Vector{Tuple{Symbol,CDorField}}
    A::Matrix{Acb}
    b::Vector{Acb}
end

"""
rels: one of `:s`, `:st`, `:su`, `:tu`, `:stu`.
      Imposes that the corresponding structure constants are
      equal, for instance rels = :st imposes that D^(s)_V = D^(t)_V.
      default: `:s` 
"""
function BootstrapMatrix(s::BootstrapSystem; rels=:s)
    unknowns, knowns = s.unknowns, s.knowns
    cache = s.fields_cache
    moduli_count = length(s.moduli_cache.s)
    if length(cache) == 3 # 3 channels
        cols_count = 0
        if rels === :stu #= [(Gs-Gt);
                             (Gs-Gu)] =#
            offset_t = 0
            offset_u = 0
            cols_count = length(unknowns.s)
            indep_channels = (:s,)
        elseif rels === :st #= [(Gs-Gt)  0 ;
                                 Gs     -Gu] =#

            offset_t = 0
            offset_u = length(unknowns.s)
            cols_count = length(unknowns.s) + length(unknowns.u)
            indep_channels = (:s,:u)
        elseif rels === :su #= [Gs        -Gt;
                                (Gs-Gu)     0] =#
            offset_t = length(unknowns.s)
            offset_u = 0
            cols_count = length(unknowns.s) + length(unknowns.t)
            indep_channels = (:s,:t)
        elseif rels === :tu #= [Gs     -Gt;
                                Gs     Gu)] =#
            offset_t = length(unknowns.s)
            offset_u = length(unknowns.s)
            cols_count = length(unknowns.s) + length(unknowns.t)
            indep_channels = (:s,:t)
        elseif rels === :s        #= [(Gs     -Gt        0
                                       Gs       0      -Gu] =#
            offset_t = length(unknowns.s)
            offset_u = length(unknowns.s) + length(unknowns.t)
            cols_count = length(unknowns.s) + length(unknowns.t) + length(unknowns.u)
            indep_channels = (:s,:t,:u)
        else
            error("rels = $rels is not supported")
        end

        A = zeros(Acb, 2*moduli_count, cols_count)
        b = zeros(Acb, 2*moduli_count)
        fields = Vector{Tuple{Symbol,CDorField}}(undef, cols_count)
        j = 0
        for chan in indep_channels
            for V in unknowns[chan]
                j += 1
                fields[j] = (chan, V)
            end
        end
        for k in 1:moduli_count
            for (j, V) in enumerate(unknowns.s)
                A[k, j]                       += cache.s[V].factors[k] * cache.s[V].block_values[k]
                A[k+moduli_count, j]          += cache.s[V].factors[k] * cache.s[V].block_values[k]
            end
            for (j, V) in enumerate(unknowns.t)
                A[k, j+offset_t]              -= cache.t[V].factors[k] * cache.t[V].block_values[k]
            end
            for (j, V) in enumerate(unknowns.u)
                A[k+moduli_count, j+offset_u] -= cache.u[V].factors[k] * cache.u[V].block_values[k]
            end

            # right-hand side of the linear system
            for (V, D) in knowns.s
                b[k] -= cache.s[V].factors[k] * D * cache.s[V].block_values[k]
                b[k + moduli_count] -= cache.s[V].factors[k] * D * cache.s[V].block_values[k]
            end
            for (V, D) in knowns.t
                b[k] += cache.t[V].factors[k] * D * cache.t[V].block_values[k]
            end
            for (V, D) in knowns.u
                b[k + moduli_count] += cache.u[V].factors[k] * D * cache.u[V].block_values[k]
            end
        end
    else
        error("BootstrapSystem with 2 channels not implemented yet")
    end
    return BootstrapMatrix(fields, A, b)
end

function solve(s::BootstrapSystem; rels=:s)
    matrix = BootstrapMatrix(s, rels=rels)
    chans = Set(matrix.fields[j][1] for j in 1:length(matrix.fields))
    str_cst_dicts = Channels(Dict(chan => Dict{CDorField,Tuple{Acb,Float64}}()
                                 for chan in chans))
    for chan in chans
        for (V, val) in s.knowns[chan]
            str_cst_dicts[chan][V] = (val, 0.0)
        end
    end

    LHS = convert(Matrix{Complex{BigFloat}}, matrix.A)
    RHS = convert(Vector{Complex{BigFloat}}, matrix.b) # convert to BigFloat to use Julia's solver

    sol1 = LHS[3:end, :] \ RHS[3:end]
    sol2 = LHS[1:(end-2), :] \ RHS[1:(end-2)]

    # compare the two solutions
    errors = @. Float64(abs((sol1 - sol2) / sol2))

    for (j, (chan, V)) in enumerate(matrix.fields)
        str_cst_dicts[chan][V] = (sol1[j], errors[j])
    end

    return Channels(Dict(chan => StructureConstants(str_cst_dicts[chan]) for chan in chans))
end

struct ChannelSpectrum
    correlation::Correlation
    blocks::Vector{Block}
    Δmax::Int
end
const ChanSpec = ChannelSpectrum

function ChannelSpectrum(co::Correlation, chan, Vs, f::Function;
                         Δmax=nothing)
    if Δmax === nothing
        Δmax = co.Δmax
    end
    blocks = Block[]
    for V in Vs
        if (real(total_dimension(V)) < Δmax)
            push!(blocks, f(V))
        end
    end
    return ChannelSpectrum(co, blocks, Δmax)
end

Base.length(s::ChanSpec) = length(s.blocks)

function BootstrapSystem(specs::Channels{ChannelSpectrum};
                         fix=nothing, moduli=nothing, extrapoints::Int=6,
                         prefactors=nothing)
    co = specs.s.correlation
    blocks = Channels(Dict(chan => specs[chan].blocks
                           for chan in keys(specs)))
    return BootstrapSystem(co, blocks;
                           fix=fix, moduli=moduli, extrapoints=extrapoints, prefactors=prefactors)
end

function solve_bootstrap(specs::Channels{ChannelSpectrum}; rels=:s)
    sys = BootstrapSystem(specs)
    solve(sys, rels=rels)
end

function to_columntable(c::SC, rmax=nothing)
    fields = sort(collect(keys(c.dict)), by = V -> real(total_dimension(V)))
    rmax !== nothing && (fields = filter(V -> V.r <= rmax, fields))
    return (
        Field = string.(fields),
        StructureConstant = [c.dict[V][1] for V in fields],
        RelativeError = [c.dict[V][2] for V in fields],
    )
end

function Base.delete!(s::SC, V)
    delete!(s.dict, V)
end

function Base.show(io::IO, c::SC; rmax = nothing, backend = :text, digits=8)
    t = columntable(to_columntable(c, rmax))
    hl_conv = Highlighter(
        (table, i, j) ->
            table[:RelativeError][i] < 2^(-precision(Acb, base = 10) / 2) &&
                j >= 2,
        crayon"green",
    )
    hl_zero = Highlighter(
        (table, i, j) ->
            table[:RelativeError][i] > 1e-2 &&
                abs(table[:StructureConstant][i]) <
                2^(-precision(Acb, base = 10) / 3) &&
                j >= 2,
        crayon"green",
    )
    hl_not_conv = Highlighter(
        (table, i, j) ->
            table[:RelativeError][i] > 1e-2 &&
                abs(table[:StructureConstant][i]) >
                2^(-precision(Acb, base = 10) / 3) &&
                j >= 2,
        crayon"red",
    )
    fmt_zero =
        (v, i, j) ->
            (
                t.RelativeError[i] > 1e-2 &&
                abs(t.StructureConstant[i]) <
                2^(-precision(Acb, base = 10) / 3) &&
                j == 2
            ) ? 0 : v
    fmt_relerr = (v, i, j) -> j == 3 && v isa Real ? format(Format("%.2g"), v) : v
    fmt_indices = (v, i, j) -> j == 1 ? tex_tuple(v) : v
    if backend == :latex
        fmt_vals = (v, i, j) -> j == 2 ? tex_complex(v) : v
        pretty_table(
            io,
            t,
            header = ["Field", "Structure constant", "Rel. err."],
            header_alignment = :l,
            formatters = (fmt_relerr, fmt_zero, fmt_indices, fmt_vals),
            backend = Val(backend),
            alignment = :l,
        )
        print(io, "\\\\")
    else
        fmt_vals = (v, i, j) -> j == 2 ? format_complex(v, digits=digits) : v
        pretty_table(
            io,
            t,
            header = ["Field", "Structure constant", "Rel. err."],
            header_alignment = :l,
            alignment = :l,
            header_crayon = crayon"blue",
            highlighters = (hl_conv, hl_zero, hl_not_conv),
            formatters = (fmt_zero, fmt_relerr, fmt_vals),
            backend = Val(backend),
            crop = :none,
        )
    end
end

function Base.show(io::IO, c::Channels{<:SC}; rmax = nothing, backend = :text)
    for chan in keys(c)
        channel_header = " Channel $chan\n"
        if backend === :latex
            print(io, channel_header)
        else
            printstyled(io, channel_header; bold = true)
        end
        show(io, c[chan], rmax = rmax, backend = backend)
    end
end
