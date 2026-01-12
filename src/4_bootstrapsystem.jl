struct StructureConstants{T<:Complex}
    constants::Dict{Field{T},T}
    errors::Dict{Field{T},Float64}
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
    channels::Tuple{Vararg{Symbol}}
    unknowns::Channels{Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
    correlation::Correlation
    positions::Vector{T}
    positions_cache::Channels{Vector{LRPosCache}} # positions at which eqs are evaluated
    spectra::Channels{ChanSpec{T}}    # channel spectra
    block_values::Channels{Dict{Field{T},Vector{T}}} # all blocks evaluated at all positions
    factors::Channels{Vector{T}} # overall position-dependent factor multiplying blocks in each channel
    str_cst::Any
    relations::Any
end

#======================================================================================
Structure constants
======================================================================================#
function StructureConstants{T}() where {T}
    constants = Dict{Field{T},T}()
    errors = Dict{Field{T},Float64}()
    return StructureConstants{T}(constants, errors)
end

function Base.getproperty(c::SC, s::Symbol)
    s === :fields && begin
        consts = getfield(c, :constants)
        return vcat([[V for V in keys(consts)]]...)
    end
    getfield(c, s)
end

Base.length(c::SC) = length(c.constants)

function fix!(cnst, field, value; error = 0)
    cnst.constants[field] = value
    cnst.errors[field] = error
    return
end

function to_dataframe(c::SC, rmax = nothing)
    df = DataFrame(
        Field = String[],
        StructureConstant = Complex[],
        RelativeError = Float64[],
    )

    fields = sort(collect(keys(c.constants)), by = V -> real(total_dimension(V)))
    rmax !== nothing && (fields = filter(V -> V.r <= rmax, fields))
    for V in fields
        push!(df, (string(V), c.constants[V], c.errors[V]))
    end
    return df
end

function Base.delete!(s::SC, V)
    delete!(s.constants, V)
    delete!(s.errors, V)
end

# Convert "(2, 1//2)" or "(1, 0)" into LaTeX math tuples
function tex_tuple(s::AbstractString)
    # Match "(a, b)" patterns where a,b may include rationals
    m = match(r"\(([^,]+),\s*([^)]+)\)", s)
    isnothing(m) && return L"$%$(s)$"  # fallback: just wrap whole string

    a = m.captures[1]
    b = m.captures[2]

    # Convert rationals inside strings, e.g. 1//2 → \frac{1}{2}
    a_tex = replace(a, r"(-?\d+)//(\d+)" => s"\\frac{\1}{\2}")
    b_tex = replace(b, r"(-?\d+)//(\d+)" => s"\\frac{\1}{\2}")

    # Wrap result in \left(...\right)
    return L"$\left(%$a_tex,%$b_tex\right)$"
end

function Base.show(io::IO, c::SC; rmax = nothing, backend = :text)
    t = columntable(to_dataframe(c, rmax))
    hl_conv = Highlighter(
        (table, i, j) ->
            table[:RelativeError][i] < 2^(-precision(BigFloat, base = 10) / 2) &&
                j >= 2,
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
    fmt_relerr = (v, i, j) -> j == 3 && v isa Real ? format(Format("%.2g"), v) : v
    fmt_indices = (v, i, j) -> j == 1 ? tex_tuple(v) : v
    fmt_vals = (v, i, j) -> j == 2 ? tex_complex(v) : v
    if backend == :latex
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
        pretty_table(
            io,
            t,
            header = ["Field", "Structure constant", "Rel. err."],
            header_alignment = :l,
            alignment = :l,
            header_crayon = crayon"blue",
            highlighters = (hl_conv, hl_zero, hl_not_conv),
            formatters = (fmt_zero, fmt_relerr),
            backend = Val(backend),
        )
    end
end

function Base.show(io::IO, c::Channels{<:SC}; rmax = nothing, backend = :text)
    for chan in CHANNELS
        channel_header = " Channel $chan\n"
        if backend === :latex
            print(io, channel_header)
        else
            printstyled(io, channel_header; bold = true)
        end
        show(io, c[chan], rmax = rmax, backend = backend)
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

function Base.:+(c1::SC{T}, c2::SC) where {T}
    consts = deepcopy(c1.constants)
    errs = deepcopy(c1.errors)
    for (k, v) in c2.constants
        if k in keys(c1.constants)
            consts[k] += v
            errs[k] += c2.errors[k]
        else
            consts[k] = v
            errs[k] = c2.errors[k]
        end
    end
    return StructureConstants{T}(consts, errs)
end

function Base.:-(c1::SC{T}, c2::SC) where {T}
    consts = deepcopy(c1.constants)
    errs = deepcopy(c1.errors)
    for (k, v) in c2.constants
        if k in keys(c1.constants)
            consts[k] -= v
            errs[k] += c2.errors[k]
        else
            consts[k] = v
            errs[k] = c2.errors[k]
        end
    end
    return StructureConstants{T}(consts, errs)
end

function find_normalised(c::SC)
    for (k, v) in c.constants
        if v ≈ 1 && c.errors[k] ≈ 0
            return k
        end
    end
    return nothing
end

function find_normalised(c::Channels)
    for chan in CHANNELS
        find = find_normalised(c[chan])
        if find  !== nothing
            return find, chan
        end
    end
    return nothing
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
    to_add = [V for V in Vs if real(total_dimension(V)) < s.Δmax]
    blocks = Vector(undef, length(to_add))
    Threads.@threads for i = eachindex(to_add)
        V = to_add[i]
        b = f(V)
        blocks[i] = b
    end
    for b in blocks
        _add_one!(s, b)
    end
    # for V in Vs
    #     if real(total_dimension(V)) < s.Δmax
    #         add!(s, f(V))
    #     end
    # end
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

function finddiag(s)
    for b in values(s.blocks)
        b.chan_field.diagonal && return b
    end
    return nothing
end

function removediag!(s)
    V = finddiag(s)
    V !== nothing && remove!(s, V)
    return V
end

function Base.filter(f, s::ChanSpec{T}) where {T}
    return ChanSpec{T}(s.corr, filter(kv -> f(kv[1]), s.blocks), s.Δmax, s.chan)
end

Base.length(s::ChanSpec) = length(s.blocks)
Base.size(s::ChanSpec) = size(s.blocks)
nb_blocks(s::ChanSpec) = sum(length(b) for b in s.blocks)

function Base.show(io::IO, s::ChanSpec)
    println(io, "channel $(s.chan), Δmax = $(s.Δmax)")
    nondiags =
        sort([(V.r, V.s) for V in fields(s) if !V.diagonal], by = x -> (x[1], x[2]))
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
    thread_results = [Dict{Field{T},Vector{T}}() for _ = 1:max_thread_id()]

    bs_vals = collect(values(bs))

    Threads.@threads for i = 1:nbs
        b = bs_vals[i]
        thread_results[Threads.threadid()][b.chan_field] = [b(x) for x in xs]
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
    rels = nothing,
) where {T}
    co = S.s.corr
    if knowns === nothing
        knowns = @channels StrCst{T}()
    end

    if rels === nothing
        channels = (:s, :t, :u)
    elseif rels === :st
        channels = (:s, :u)
        knowns.t = knowns.s
    elseif rels === :su
        channels = (:s, :t)
        knowns.u = knowns.s
    elseif rels === :tu
        channels = (:s, :t)
        knowns.u = knowns.t
    elseif rels === :stu
        channels = (:s,)
        knowns.t = knowns.s
        knowns.u = knowns.s
    else
        error("rels has to be one of :stu, :st, :su, :tu, nothing")
    end
    nb_channels = length(channels)
    nb_unknowns = [length(S[chan]) for chan in channels]
    nb_unknowns .-= length.([knowns[chan] for chan in channels])
    nb_lines = 2 * ((sum(nb_unknowns) + extrapoints) ÷ 2)
    if pos !== nothing
        nb_positions = length(pos)
        @assert nb_positions > nb_lines ÷ (nb_channels - 1) "
           You must give enough positions for the system to be overdetermined
       "
    else
        nb_positions = nb_channels == 1 ? nb_lines ÷ 2 : nb_lines ÷ (nb_channels - 1)
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
    factors = @channels [factor(co, chan, p) for p in pos]

    return BootstrapSystem{T}(
        channels,
        unknowns,
        LHS,
        RHS,
        co,
        pos,
        pos_cache,
        S,
        blocks,
        factors,
        knowns,
        rels,
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

function computeLHScolumn(b::BootstrapSystem{T}, chan, V) where {T}
    rels = b.relations
    zer = zeros(T, length(b.positions))
    if rels === nothing
        v = b.factors[chan] .* b.block_values[chan][V]
        chan === :s && return vcat(v, v)
        chan === :t && return vcat(-v, zer)
        chan === :u && return vcat(zer, -v)
    elseif rels == :st  #= [(Gs-Gt)  0;
                                  Gs     -Gu] =#
        if chan == :s
            vs_s = b.factors[:s] .* b.block_values[:s][V]
            vs_t = b.factors[:t] .* b.block_values[:t][V]
            return vcat(vs_s .- vs_t, vs_s)
        else
            vs_u = b.factors[:u] .* b.block_values[:u][V]
            return vcat(zer, .-vs_u)
        end
    elseif rels == :su  #= [Gs     -Gt;
                                  (Gs-Gu)   0] =#
        if chan == :s
            vs_s = b.factors[:s] .* b.block_values[:s][V]
            vs_u = b.factors[:u] .* b.block_values[:u][V]
            return vcat(vs_s, vs_s .- vs_u)
        else
            vs_t = b.factors[:t] .* b.block_values[:t][V]
            return vcat(.-vs_t, zer)
        end
    elseif rels == :tu   #= [Gs     -Gt;
                                  0   (Gt-Gu)] =#
        if chan == :s
            vs_s = b.factors[:s] .* b.block_values[:s][V]
            return vcat(vs_s, zer)
        else
            vs_t = b.factors[:t] .* b.block_values[:t][V]
            vs_u = b.factors[:u] .* b.block_values[:u][V]
            return vcat(.-vs_t, vs_t .- vs_u)
        end
    elseif rels == :stu #= [(Gs-Gt);
                                   (Gs-Gu)] =#
        vs = @channels b.factors[chan] .* b.block_values[chan][V]
        return vcat(vs.s .- vs.t, vs.s .- vs.u)
    end
end

function computeRHS(b::BootstrapSystem{T}, knowns) where {T}
    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    nb_positions = length(b.positions)
    sumknowns(knowns, chan, i) = sum(
        b.str_cst[chan].constants[V] *
        b.factors[chan][i] *
        b.block_values[chan][V][i] for V in knowns[chan];
        init = zero(T),
    )
    return vcat(
        [
            sumknowns(knowns, :t, i) - sumknowns(knowns, :s, i) for
            i in eachindex(b.positions)
        ],
        [
            sumknowns(knowns, :u, i) - sumknowns(knowns, :s, i) for
            i in eachindex(b.positions)
        ],
    )
end

function compute_linear_system!(b::BootstrapSystem{T}) where {T}
    known_consts = b.str_cst

    unknowns = @channels sort(
        [
            V for V in fields(b.spectra[chan]) if
            !(V in keys(known_consts[chan].constants))
        ],
        by = V -> real(total_dimension(V)),
    )
    for chan in CHANNELS
        if !(chan in b.channels)
            empty!(unknowns[chan])
        end
    end
    # if no consts are known, fix a normalisation
    if all([isempty(known_consts[chan].constants) for chan in b.channels])
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
        fix!(known_consts[chan], normalised_field, 1)
        b.str_cst = known_consts
        idx = findfirst(==(normalised_field), unknowns[chan])
        deleteat!(unknowns[chan], idx)
    end

    knowns = @channels [
        V for V in fields(b.spectra[chan]) if V in keys(b.str_cst[chan].constants)
    ]

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(CHANNELS) - 1)
    nb_unknowns = sum(
        length(b.spectra[chan]) - length(known_consts[chan]) for chan in b.channels
    )

    # form the matrix of the linear system
    # keep a record of the columns where we put each field.
    LHS = Matrix{T}(undef, nb_lines, nb_unknowns)
    col_idx = 0
    for chan in b.channels
        for V in unknowns[chan]
            col_idx += 1
            LHS[:, col_idx] = computeLHScolumn(b, chan, V)
        end
    end

    RHS = computeRHS(b, knowns)
    b.unknowns = unknowns
    b.LHS = LHS
    b.RHS = RHS
    return
end

function findcol(s::BootstrapSystem, V::Field, chan)
    col = 1
    uk = s.unknowns
    chan === :t && (col += length(uk.s)) # t blocks start at this col
    chan === :u && (col += length(uk.s) + length(uk.t))
    return col - 1 + findfirst(isequal(V), uk[chan])
end

function replacediagLHS!(s::BootstrapSystem, new::Block)
    # locate the column of the block to remove.
    chan = new.chan
    old = removediag!(s.spectra[chan])
    remove!(s.spectra[chan], old.chan_field)
    add!(s.spectra[chan], new)
    newV = new.chan_field
    col = findcol(s, old.chan_field, chan)
    s.unknowns[chan][findfirst(==(old.chan_field), s.unknowns[chan])] = newV
    delete!(s.block_values[chan], old.chan_field)
    s.block_values[chan][newV] = [new(x) for x in s.positions_cache[chan]]
    s.LHS[:, col] = computeLHScolumn(s, chan, new.chan_field)
end

function replacediagRHS!(sys::BootstrapSystem, new::Block, str_cst = 1)
    chan = new.chan
    old = removediag!(sys.spectra[chan])
    remove!(sys.spectra[chan], old.chan_field)
    add!(sys.spectra[chan], new)
    newV = new.chan_field
    oldV = old.chan_field
    sys.block_values[chan][newV] = [new(x) for x in sys.positions_cache[chan]]
    to_add =
        sys.str_cst[chan].constants[oldV] .* sys.block_values[chan][oldV] .-
        str_cst .* sys.block_values[chan][newV]
    if chan === :s
        sys.RHS[1:length(sys.positions)] .+= to_add
        sys.RHS[length(sys.positions)+1:end] .+= to_add
    elseif chan === :t
        sys.RHS[1:length(sys.positions)] .-= to_add
    elseif chan === :u
        sys.RHS[length(sys.positions)+1:end] .-= to_add
    end
    delete!(sys.block_values[chan], oldV)
    delete!(sys.str_cst[chan], oldV)
    fix!(sys.str_cst[chan], newV, str_cst)
end

function replacediag!(sys::BootstrapSystem, new::Block, str_cst = 1)
    # if sys has a diagonal field with unknown str_cst in this channel, change LHS
    # else change RHS
    for V in sys.unknowns[new.chan]
        if V.diagonal
            replacediagLHS!(sys, new)
            return
        end
    end
    replacediagRHS!(sys, new, str_cst)
    return
end

function solve!(s::BootstrapSystem)
    # solve for two different sets of positions
    sol1 = s.LHS[3:end, :] \ s.RHS[3:end]
    sol2 = s.LHS[1:(end-2), :] \ s.RHS[1:(end-2)]

    # compare the two solutions
    errors = @. Float64(abs((sol1 - sol2) / sol2))

    for chan in s.channels
        for V in s.unknowns[chan]
            index = findcol(s, s.spectra[chan].blocks[V].chan_field, chan)
            s.str_cst[chan].constants[V] = sol1[index]
            s.str_cst[chan].errors[V] = errors[index]
        end
    end

    return s.str_cst
end

function eval_blocks_compute_system(specs::Channels; fix = nothing, rels = nothing)
    sys = BootstrapSystem(specs, rels = rels)
    evaluate_blocks!(sys)
    if fix !== nothing
        for (chan, V, cst) in fix
            fix!(sys.str_cst[chan], V, cst)
        end
    end
    compute_linear_system!(sys)
    return sys
end

function solve_bootstrap(specs::Channels; fix = nothing, rels = nothing)
    sys = eval_blocks_compute_system(specs, fix = fix, rels = rels)
    solve!(sys)
    return sys
end

# expects diag_blocks as a struct Channels(collection of diagonal blocks)
# if several channels have a non-empty list, the lists must be identical.
function solve_bootstrap(specs::Channels, diag_blocks; fix = nothing, rels = nothing)
    sys = eval_blocks_compute_system(specs, fix = fix, rels = rels)
    res = []
    chans = []
    if all([isempty(diag_blocks[chan]) for chan in CHANNELS])
        error(
            "You need to provide diagonal blocks in at least one of the channels to use this method",
        )
    end
    bs = []
    for chan in CHANNELS
        if !isempty(diag_blocks[chan])
            bs = diag_blocks[chan]
            push!(chans, chan)
        end
    end
    for i in eachindex(diag_blocks[chans[1]])
        for chan in chans
            replacediag!(sys, diag_blocks[chan][i])
        end
        solve!(sys)
        push!(res, deepcopy(sys))
    end
    return res
end

function clear!(b::BootstrapSystem{T}) where {T}
    b.str_cst = @channels StrCst(b.spectra[chan])
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
