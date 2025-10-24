struct StructureConstants{T<:Complex}
    constants::Channels{Dict{Field{T},T}}
    errors::Channels{Dict{Field{T},Float32}}
end
const StrCst = StructureConstants
const SC = StructureConstants

struct ChannelSpectrum{T}
    corr::Corr
    blocks::Dict{Field{T},Block{T}}
    Δmax::Int
    chan::Symbol
end
const ChanSpec = ChannelSpectrum

mutable struct BootstrapSystem{T}
    unknowns::Channels{Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
    correlation::Correlation
    positions::Vector{T}
    positions_cache::Channels{Vector{LRPosCache}} # positions at which eqs are evaluated
    spectra::Channels{ChanSpec{T}}    # channel spectra
    block_values::Channels{Dict{Field{T},Vector{T}}} # all blocks evaluated at all positions
    str_cst::StrCst{T}
end

#======================================================================================
Structure constants
======================================================================================#
function StructureConstants{T}() where {T}
    constants = @channels Dict{Field{T},T}()
    errors = @channels Dict{Field{T},Float32}()
    return StructureConstants{T}(constants, errors)
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

Base.length(c::SC) = sum(length(c.constants[chan]) for chan in keys(c.constants))

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
        t = columntable(select(table, Not(:Channel)))
        hl_conv = Highlighter(
            (table, i, j) ->
                table[:RelativeError][i] < 2^(-precision(BigFloat, base = 10) / 2) && j >= 2,
            crayon"green",
        )
        hl_zero = Highlighter(
            (table, i, j) ->
                table[:RelativeError][i] > 1e-2 &&
                abs(table[:StructureConstant][i]) <
                2^(-precision(BigFloat, base = 10) / 3) &&
                j >= 2,
            crayon"green",
        )
        hl_not_conv = Highlighter(
            (table, i, j) ->
                table[:RelativeError][i] > 1e-2 &&
                abs(table[:StructureConstant][i]) >
                2^(-precision(BigFloat, base = 10) / 3) &&
                j >= 2,
            crayon"red",
        )
        fmt_zero =
            (v, i, j) ->
                (
                    t.RelativeError[i] > 1e-2 &&
                    abs(t.StructureConstant[i]) <
                    2^(-precision(BigFloat, base = 10) / 3) &&
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

#======================================================================================
Channel spectra
======================================================================================#
fields(s::ChanSpec) = keys(s.blocks)

function _add_one!(s::ChanSpec, b)
    if !(b.chan_field in fields(s)) &&
       real(total_dimension(b.chan_field)) < real(s.Δmax)
        (s.blocks)[b.chan_field] = b
        return b
    end
end

_add!(s::ChanSpec, b::Block) = _add_one!(s, b)
_add!(s::ChanSpec, blocks::Vector) = foreach(b -> _add_one!(s, b), blocks)
add!(s, bs) = _add!(s, bs)

function add(s, blocks)
    s2 = deepcopy(s)
    add!(s2, blocks)
    return s2
end

function hasdiagonals(s::ChanSpec)
    for V in keys(s.blocks)
        V.diagonal && return true
    end
    return false
end

function ChannelSpectrum(co::Co{T}, chan, Vs::Vector{<:Field}, f::Function) where {T}
    s = ChannelSpectrum{T}(co, Dict{Field{T},Block{T}}(), co.Δmax, chan)
    for V in Vs
        if real(total_dimension(V)) < s.Δmax
            add!(s, f(V))
        end
    end
    return s
end

ChannelSpectrum(co::Co, Vs::Vector{<:Field}, f::Function) =
    ChannelSpectrum(co, :τ, Vs, f)

function _remove_one!(s::ChanSpec, V)
    delete!(s.blocks, V)
end

_remove!(s::ChanSpec, V) = _remove_one!(s, V)
_remove!(s::ChanSpec, fields::Vector) = foreach(V -> _remove_one!(s, V), fields)
remove!(s, V) = _remove!(s, V)
remove(s, Vs) = remove!(deepcopy(s), Vs)

function Base.filter(f, s::ChanSpec{T}) where {T}
    return ChanSpec{T}(s.corr, filter(kv -> f(kv[1]), s.blocks), s.Δmax, s.chan)
end

Base.length(s::ChanSpec) = length(s.blocks)
Base.size(s::ChanSpec) = size(s.blocks)
nb_blocks(s::ChanSpec) = sum(length(b) for b in s.blocks)

function Base.show(io::IO, s::ChanSpec)
    println(io, "channel $(s.chan), Δmax = $(s.Δmax)")
    nondiags =
        sort([indices(V) for V in fields(s) if !V.diagonal], by = x -> (x[1], x[2]))
    diags = sort(
        [V for V in fields(s) if V.diagonal && !V.degenerate],
        by = V -> real(total_dimension(V)),
    )
    degs = sort(
        [V for V in fields(s) if V.degenerate],
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

function evaluate_blocks(S::ChanSpec, xs)
    bs = S.blocks
    T = typeof(bs).parameters[1].parameters[1]
    nbs = length(bs)

    # Preallocate thread-local results
    thread_results = [Dict{Field{T},Vector{T}}() for _ = 1:Threads.nthreads()]

    bs_vals = collect(values(bs))

    Threads.@threads for i = 1:nbs
        b = bs_vals[i]
        tid = Threads.threadid()
        thread_results[tid][b.chan_field] = [b(x) for x in xs]
    end

    # Merge all thread-local dictionaries into a single result
    results = Dict{Field{T},Vector{T}}()
    for d in thread_results
        merge!(results, d)
    end

    return results
end

#======================================================================================
Bootstrap System
======================================================================================#
function evaluate_blocks!(sys::BootstrapSystem)
    sys.block_values =
        @channels evaluate_blocks(sys.spectra[chan], sys.positions_cache[chan])
    return sys.block_values
end

function new_random_point!(random_points, N; square = true, cond = p -> true)
    xmin, xmax, sep = square ? (0.1, 0.5, 0.2) : (-0.4, 1.4, 0.4)
    while true
        z = xmin + (xmax - xmin) * rand() + (1 + 4 * rand()) * im / 10
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
    Random.seed!(1234)
    res = []
    for _ = 1:N
        new_random_point!(res, N, square = square, cond = cond)
    end
    return transfo !== nothing ? transfo.(res) : res

    # Restore RNG state
    Random._GLOBAL_RNG[] = rng
end

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

channel_position(x, _::Correlation4, chan) = crossratio(chan, x)
channel_position(x, _::Correlation1, chan) = modular_param(chan, x)

function choose_block_eval_points(npositions, _::Correlation1{T}) where {T}
    Vector{T}(random_points(npositions, transfo = x -> τfromx(x)))
end
function choose_block_eval_points(npositions, _::Correlation4{T}) where {T}
    Vector{T}(random_points(npositions, transfo = nothing))
end

function BootstrapSystem(
    S::Channels{ChanSpec{T}};
    knowns = nothing,
    extrapoints::Int = 6,
    pos = nothing,
) where {T}
    co = S.s.corr
    if knowns === nothing
        knowns = StrCst{T}()
    end
    nb_unknowns = [length(S[chan]) for chan in CHANNELS]
    nb_unknowns .-= length.([knowns[chan] for chan in CHANNELS])
    nb_lines = 2 * ((sum(nb_unknowns) + extrapoints) ÷ 2)
    if pos !== nothing
        nb_positions = length(pos)
        @assert nb_positions > nb_lines ÷ (NB_CHANNELS - 1) "
           You must give enough positions for the system to be overdetermined
       "
    else
        nb_positions = nb_lines ÷ (NB_CHANNELS - 1)
        pos = choose_block_eval_points(nb_positions, co)
    end

    # compute position cache
    pos_chan = @channels [channel_position(x, co, chan) for x in pos]
    pos_cache = @channels Vector{LRPosCache}(undef, length(pos_chan[chan])) # result container
    Threads.@threads for i = 1:NB_CHANNELS
        chan = CHANNELS[i]
        pos_cache[chan] =
            [LRPosCache(x, S[CHANNELS[i]].corr, chan) for x in pos_chan[chan]]
    end

    # blocks = evaluate_blocks(S, pos_cache)
    blocks = @channels Dict{Field{T},Vector{T}}()
    unknowns = @channels Vector{Field{T}}()
    LHS = Matrix{T}(undef, 0, 0)
    RHS = Vector{T}()

    return BootstrapSystem{T}(
        unknowns,
        LHS,
        RHS,
        co,
        pos,
        pos_cache,
        S,
        blocks,
        knowns,
    )
end

factor(_::Correlation4{T}, chan, x) where {T} = one(T)
function factor(co::Correlation1{T}, chan, τ) where {T}
    chan === :s && return one(T)
    Δ₁, bΔ₁ = (co.fields[1].dims[lr].Δ for lr in (:left, :right))
    chan === :t && return τ^-Δ₁ * conj(τ)^-bΔ₁
    chan === :u && return (τ - 1)^-Δ₁ * conj(τ - 1)^-bΔ₁
end

function hasdiagonal(b::BootstrapSystem)
    for chan in CHANNELS
        for V in keys(b.spectra[chan].blocks)
            V.diagonal && return true, chan, V
        end
    end
    return (false, :s, nothing)
end

function compute_linear_system!(b::BootstrapSystem{T}) where {T}
    known_consts = b.str_cst
    facts = @channels [factor(b.correlation, chan, pos) for pos in b.positions]

    unknowns = @channels sort(
        [V for V in fields(b.spectra[chan]) if !(V in keys(known_consts[chan]))],
        by = V -> real(total_dimension(V)),
    )

    # if no consts are known, fix a normalisation
    if all([isempty(known_consts[chan]) for chan in CHANNELS])
        hasdiag, chan, diag = hasdiagonal(b)
        if hasdiag
            normalised_field = diag
        else
            normalised_field = reduce(
                (V1, V2) ->
                    real(total_dimension(V1)) < real(total_dimension(V2)) ? V1 : V2,
                unknowns[chan],
            )
        end
        fix!(known_consts, chan, normalised_field, 1)
        b.str_cst = known_consts
        idx = findfirst(==(normalised_field), unknowns[chan])
        deleteat!(unknowns[chan], idx)
    end

    knowns =
        @channels [V for V in fields(b.spectra[chan]) if V in keys(b.str_cst[chan])]

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(CHANNELS) - 1)
    nb_unknowns = [length(b.spectra[chan]) for chan in CHANNELS]
    nb_unknowns .-= length.([known_consts[chan] for chan in CHANNELS])

    # matrix of equations Σ_unknowns (chan1 - chan2) = Σ_knowns (chan2 - chan1)
    # and Σ_unknowns (chan1 - chan3) = Σ_knowns (chan3 - chan1)
    zer = zeros(T, nb_positions)

    # Form the matrix
    # [ G^s -G^t  0  ;
    #   G^s  0   -G^u ]
    LHS = Matrix{T}(undef, nb_lines, sum(nb_unknowns))
    col_idx = 1
    for V in unknowns.s
        LHS[1:nb_positions, col_idx] = facts.s .* b.block_values.s[V]
        LHS[(nb_positions+1):end, col_idx] = facts.s .* b.block_values.s[V]
        col_idx += 1
    end
    for V in unknowns.t
        LHS[1:nb_positions, col_idx] = .-facts.t .* b.block_values.t[V]
        LHS[(nb_positions+1):end, col_idx] = zer
        col_idx += 1
    end
    for V in unknowns.u
        LHS[1:nb_positions, col_idx] = zer
        LHS[(nb_positions+1):end, col_idx] = .-facts.u .* b.block_values.u[V]
        col_idx += 1
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = Vector{T}(undef, 2nb_positions)
    for i = 1:nb_positions
        RHS[i] = sum(b.block_values.t[V][i] for V in knowns.t; init = zero(T))
        RHS[i] -= sum(b.block_values.s[V][i] for V in knowns.s; init = zero(T))
        RHS[i+nb_positions] =
            sum(b.block_values.u[V][i] for V in knowns.u; init = zero(T))
        RHS[i+nb_positions] -=
            sum(b.block_values.s[V][i] for V in knowns.s; init = zero(T))
    end

    b.unknowns = unknowns
    b.LHS = LHS
    b.RHS = RHS
end

function solve!(s::BootstrapSystem; precision_factor = 1)
    # solve for two different sets of positions
    if precision_factor > 1
        s.LHS = Matrix{Complex{BigFloat}}(s.LHS)
        s.RHS = Vector{Complex{BigFloat}}(s.RHS)
    end
    sol1, sol2 = setprecision(BigFloat, precision_factor * precision(BigFloat)) do
        (s.LHS[3:end, :] \ s.RHS[3:end], s.LHS[1:(end-2), :] \ s.RHS[1:(end-2)])
    end

    # compare the two solutions
    errors = @. Float32(abs((sol1 - sol2) / sol2))

    # back to dictionary format
    nb_unknowns = vcat([0], [length(s.unknowns[chan]) for chan in CHANNELS])
    for (i, chan) in enumerate(CHANNELS)
        for (j, V) in enumerate(s.unknowns[chan])
            index = sum(nb_unknowns[k] for k = 1:(i-1); init = 0) + j
            s.str_cst[chan][V] = sol1[index]
            s.str_cst.errors[chan][V] = errors[index]
        end
    end

    return s.str_cst
end

function solve_bootstrap(specs::Channels)
    sys = BootstrapSystem(specs)
    evaluate_blocks!(sys)
    compute_linear_system!(sys)
    solve!(sys)
    return sys
end

function clear!(b::BootstrapSystem{T}) where {T}
    b.str_cst = StrCst(b.spectra)
    # b.matrix = BootstrapMatrix{T}([chan for chan in keys(b.spectra)])
end

function evaluate_correlation(sys, chan)
    z = random_points(1)[1]
    co = sys.correlation
    consts = sys.str_cst
    cache = LRPosCache(channel_position(z, co, chan), co, chan)
    block_values =
        Dict(V => b(cache) for (V, b) in consts[chan] if consts.errors[chan] < 1e-4)
    return sum(consts[chan][V] * block_values[V] for V in keys(block_values))
end

function Base.show(io::IO, b::BootstrapSystem)
    println(io, "BootstrapSystem{$(typeof(b).parameters[1])}")
    println(io, b.correlation)
    println(io, "Number of positions: $(length(b.positions))")
    println(io, "Matrix size: $(size(b.LHS))")
end
