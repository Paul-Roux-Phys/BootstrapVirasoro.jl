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

# function BootstrapMatrix{T}() where {T}
#     BootstrapMatrix{T}(
#         Channels(Vector{Field{T}}())Matrix{T}(undef, 0, 0),
#         Vector{T}(undef, 0),
#     )
# end

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

evaluate_blocks(S::Channels{ChanSpec}, xs) =
    @channels evaluate_blocks(S[chan], xs[chan])

function evaluate_blocks!(sys::BootstrapSystem)
    sys.block_values = evaluate_blocks(sys.spectra, sys.positions_cache)
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
    res = []
    for _ = 1:N
        new_random_point!(res, N, square = square, cond = cond)
    end
    return transfo !== nothing ? transfo.(res) : res
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
    chans = (:s, :t, :u)
    co = S[1].corr
    if knowns === nothing
        knowns = StrCst{T}()
    end
    nb_unknowns = [length(s) for s in S]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])
    nb_lines = 2 * ((sum(nb_unknowns) + extrapoints) ÷ 2)
    if pos !== nothing
        nb_positions = length(pos)
        @assert nb_positions > nb_lines ÷ (NB_CHANNELS-1) "
            You must give enough positions for the system to be overdetermined
        "
    else
        nb_positions = nb_lines ÷ (NB_CHANNELS-1)
        pos = choose_block_eval_points(nb_positions, co)
    end

    # compute position cache
    pos_chan = @channels [channel_position(x, co, chan) for x in pos]
    pos_cache = @channels Vector{LRPosCache{T}}(undef, length(pos_chan[chan])) # result container
    Threads.@threads for i in 1:NB_CHANNELS
        chan = chans[i]
        pos_cache[chan] = [LRPosCache(x, S[i].corr, chan) for x in pos_chan[chan]]
    end

    # blocks = evaluate_blocks(S, pos_cache)
    blocks = (; (chan => Dict{Field{T},Vector{T}}() for chan in chans)...)
    unknowns =  @channels Vector{Field{T}}()
    LHS = Matrix{T}(undef, 0, 0)
    RHS = Vector{T}()

    return BootstrapSystem{T}(unknowns, LHS, RHS, co, pos, pos_cache, S, blocks, knowns)
end

function add!(sys::BootstrapSystem, chan, V)
    new_block = add!(sys.spectra[chan], V)
    if isnothing(new_block)
        error("trying to add already present block")
        return
    end
    # add a new position
    new_random_point!(sys.positions, length(sys.positions))
    new_pos = sys.positions[end]
    for chan in CHANNELS
        new_pos_cache = LRPosCache(new_pos, sys.spectra[chan].corr, chan)
        push!(sys.positions_cache[chan], new_pos_cache)
    end
    for chan in CHANNELS
        # evaluate previous blocks at the new position
        new_pos = sys.positions_cache[chan][end]
        for v in sys.block_values[chan]
            push!(v, v(new_pos))
        end
    end
    # evaluate the new block at all positions
    sys.block_values[chan][V] = [new_block(x) for x in sys.positions_cache[chan]]
end

function remove!(b::BootstrapSystem, chan, V)
    remove!(b.spectra[chan], V)
    pop!(b.positions)
    for chan in CHANNELS
        pop!(b.positions_cache[chan])
    end
    # remove one position for all blocks
    for chan in CHANNELS
        for v in b.block_values[chan]
            pop!(v)
        end
    end
    # remove the block
    delete!(b.block_values[chan], V)
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
    facts = Dict(
        chan => [factor(b.correlation, chan, pos) for pos in b.positions] for
        chan in CHANNELS
    )

    unknowns = (;
        (
            chan => sort(
                [
                    V for V in fields(b.spectra[chan]) if
                    !(V in keys(known_consts[chan]))
                ],
                by = V -> real(total_dimension(V)),
            ) for chan in CHANNELS
        )...
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

    knowns = (;
        (
            chan => [
                V for V in fields(b.spectra[chan]) if
                V in keys(b.str_cst[chan])
            ] for chan in CHANNELS
        )...
    )

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(CHANNELS) - 1)
    nb_unknowns = [length(s) for s in b.spectra]
    nb_unknowns .-= length.([known_consts[chan] for chan in CHANNELS])

    # matrix of equations Σ_unknowns (chan1 - chan2) = Σ_knowns (chan2 - chan1)
    # and Σ_unknowns (chan1 - chan3) = Σ_knowns (chan3 - chan1)
    zer = zeros(T, nb_positions)

    # Form the matrix
    # [ G^s -G^t  0  ;
    #   G^s  0   -G^u ]
    LHS = Matrix{T}(undef, nb_lines, sum(nb_unknowns))
    col_idx = 1
    for V in unknowns[:s]
        LHS[1:nb_positions, col_idx] = facts[:s] .* b.block_values[:s][V]
        LHS[(nb_positions+1):end, col_idx] = facts[:s] .* b.block_values[:s][V]
        col_idx += 1
    end
    for V in unknowns[:t]
        LHS[1:nb_positions, col_idx] = .-facts[:t] .* b.block_values[:t][V]
        LHS[(nb_positions+1):end, col_idx] = zer
        col_idx += 1
    end
    for V in unknowns[:u]
        LHS[1:nb_positions, col_idx] = zer
        LHS[(nb_positions+1):end, col_idx] = .-facts[:u] .* b.block_values[:u][V]
        col_idx += 1
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = Vector{T}(undef, 2nb_positions)
    for i = 1:nb_positions
        RHS[i] = sum(b.block_values[:t][V][i] for V in knowns[:t]; init = zero(T))
        RHS[i] -= sum(b.block_values[:s][V][i] for V in knowns[:s]; init = zero(T))
        RHS[i+nb_positions] =
            sum(b.block_values[:u][V][i] for V in knowns[:u]; init = zero(T))
        RHS[i+nb_positions] -=
            sum(b.block_values[:s][V][i] for V in knowns[:s]; init = zero(T))
    end

    b.unknowns = unknowns
    b.LHS = LHS
    b.RHS = RHS
end

function solve!(s::BootstrapSystem; precision_factor = 1)
    # solve for two different sets of positions
    LHS, RHS = s.LHS, s.RHS
    if precision_factor > 1
        LHS = Matrix{Complex{BigFloat}}(LHS)
        RHS = Vector{Complex{BigFloat}}(RHS)
    end
    sol1, sol2 = setprecision(BigFloat, precision_factor * precision(BigFloat)) do
        (
            s.LHS[3:end, :] \ s.RHS[3:end],
            s.LHS[1:(end-2), :] \ s.RHS[1:(end-2)],
        )
    end

    # compare the two solutions
    errors = @. Float32(abs((sol1 - sol2) / sol2))

    # back to dictionary format
    chans = keys(s.spectra)
    nb_unknowns = vcat([0], [length(s.unknowns[chan]) for chan in chans])
    for (i, chan) in enumerate(chans)
        for (j, V) in enumerate(s.unknowns[chan])
            index = sum(nb_unknowns[k] for k = 1:(i-1); init = 0) + j
            s.str_cst[chan][V] = sol1[index]
            s.str_cst.errors[chan][V] = errors[index]
        end
    end

    return s.str_cst
end

function clear!(b::BootstrapSystem{T}) where {T}
    b.str_cst = StrCst(b.spectra)
    # b.matrix = BootstrapMatrix{T}([chan for chan in keys(b.spectra)])
end

compute_reference!(sys::BootstrapSystem) =
    compute_reference!(sys.str_cst, sys.spectra)

function evaluate_correlation(sys, chan)
    z = random_points(1)[1]
    co = sys.correlation
    consts = sys.str_cst
    cache = LRPosCache(channel_position(z, co, chan), co, chan)
    block_values =
        Dict(V => b(cache) for (V, b) in consts[chan] if consts.errors[chan] < 1e-4)
    return sum(consts[chan][V] * block_values[V] for V in keys(block_values))
end
