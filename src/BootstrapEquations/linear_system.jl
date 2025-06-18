struct BootstrapMatrix{T}
    unknowns::Channels{Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
end

"""
TODO
"""
mutable struct BootstrapSystem{T,U<:ChannelSpectrum{T}}
    correlation::Correlation
    positions::Vector{T}
    positions_cache::Channels{Vector{LRPositionCache{T}}}        # positions at which eqs are evaluated
    spectra::Channels{U}    # channel spectra
    block_values::Channels{Dict{Field{T}, Vector{T}}} # all blocks evaluated at all positions
    matrix::BootstrapMatrix{T}  # matrix of equations
    consts::StructureConstants{T}
end

function BootstrapMatrix{T}() where {T}
    BootstrapMatrix{T}(
        Channels{Vector{Field{T}}}(
            Tuple(Vector{Field{T}}() for chan in (:s, :t, :u)),
        ),
        Matrix{T}(undef, 0, 0),
        Vector{T}(undef, 0),
    )
end

function evaluate_blocks(bs::Dict, xs)
    T = typeof(bs).parameters[1].parameters[1]
    nbs = length(bs)

    # Preallocate thread-local results
    thread_results = [Dict{Field{T}, Vector{T}}() for _ in 1:Threads.nthreads()]

    bs_vals = collect(values(bs))

    Threads.@threads for i in 1:nbs
        b = bs_vals[i]
        tid = Threads.threadid()
        thread_results[tid][b.channel_field] = [b(x) for x in xs]
    end

    # Merge all thread-local dictionaries into a single result
    results = Dict{Field{T}, Vector{T}}()
    for d in thread_results
        merge!(results, d)
    end

    return results
end

function evaluate_blocks(S::Channels{U}, xs) where {T,U<:ChannelSpectrum{T}}
    return Channels{Dict{Field{T}, Vector{T}}}(Tuple(
        evaluate_blocks(S[chan].blocks, xs[chan]) for chan in keys(S)
    ))
end

function evaluate_blocks!(sys::BootstrapSystem)
    sys.block_values = evaluate_blocks(sys.spectra, sys.positions_cache)
    return sys.block_values
end

function new_random_point!(random_points, N; transfo=missing, square=true)
    xmin, xmax, sep = square ? (0.1, 0.5, 0.2) : (-0.4, 1.4, 0.4)
    while true
        z = xmin + (xmax - xmin) * rand() + (1 + 4 * rand()) * im / 10
        if minimum(append!(abs.(z .- random_points), 1)) > sep / sqrt(N)
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
function random_points(N; transfo=missing, square=true)
    res = []
    for _ = 1:N
        new_random_point!(res, N, transfo=transfo, square=square)
    end
    return transfo !== missing ? transfo.(res) : res
end

function BootstrapSystem(
    S::Channels{U};
    knowns=nothing,
    extrapoints::Int=6,
    pos=nothing,
) where {T,U<:ChannelSpectrum{T}}
    # if consts is empty, normalise the field with smallest indices
    chans = keys(S)
    co = S[1].corr
    if knowns === nothing
        knowns = StructureConstants{T}()
    end
    nb_unknowns = [length(s) for s in S]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])
    nb_lines = 2 * ((sum(nb_unknowns) + extrapoints) ÷ 2)
    if pos !== nothing
        nb_positions = length(pos)
    else
        nb_positions = nb_lines ÷ (length(S) - 1)
        pos = Vector{T}(random_points(nb_positions))
    end
    pos_chan = Channels{Vector{T}}(
        Tuple([channel_position(x, co, chan) for x in pos] for chan in keys(S)),
    )
    # compute position cache
    pos_cache_data = Vector{Vector{LRPositionCache{T}}}(undef, length(S)) # result container
    Threads.@threads for i = 1:3
        chan = chans[i]
        pos_cache_data[i] = [LRPositionCache(x, S[i].corr, chan) for x in pos_chan[i]]
    end
    pos_cache = Channels{Vector{LRPositionCache{T}}}(Tuple(pos_cache_data))
    # blocks = evaluate_blocks(S, pos_cache)
    blocks = Channels{Dict{Field{T}, Vector{T}}}(Tuple(Dict{Field{T}, Vector{T}}() for chan in chans))
    matrix = BootstrapMatrix{T}()
    return BootstrapSystem{T,U}(co, pos, pos_cache, S, blocks, matrix, knowns)
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
    for chan in (:s, :t, :u)
        new_pos_cache = LRPositionCache(new_pos, sys.spectra[chan].corr, chan)
        push!(sys.positions_cache[chan], new_pos_cache)
    end
    for chan in (:s, :t, :u)
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
    for chan in (:s, :t, :u)
        pop!(b.positions_cache[chan])
    end
    # remove one position for all blocks
    for chan in (:s, :t, :u)
        for v in b.block_values[chan]
            pop!(v)
        end
    end
    # remove the block
    delete!(b.block_values[chan], V)
end

factor(_::Correlation{T, U}, chan, x) where {T, U<:FourPoints} = one(T)
function factor(co::Correlation{T, U}, chan, τ) where {T, U<:OnePoint}
    Δ₁, bΔ₁ = (co.fields[1].dims[lr].Δ for lr in (:left, :right))
    chan === :t && return (-1)^(Δ₁) * τ^(-Δ₁) * τ^(-bΔ₁)
    return one(T)
end

function compute_linear_system!(b::BootstrapSystem{T}) where {T}
    chans = (:s, :t, :u)
    known_consts = b.consts
    facts = Dict(
        chan => [factor(b.correlation, chan, pos) for pos in b.positions]
        for chan in chans
    )

    unknowns = Channels{Vector{Field{T}}}(Tuple(sort([
        V for V in fields(b.spectra[chan])
        if !(V in keys(known_consts[chan]))
            ], by=V->real(total_dimension(V))) for chan in chans
    ))

    # if no consts are known, fix a normalisation: first str cst in the s channel = 1
    if all([isempty(known_consts[chan]) for chan in chans])
        smallest_dim = reduce(
            (V1, V2) -> real(total_dimension(V1)) < real(total_dimension(V2)) ? V1 : V2, unknowns[:s]
        )
        fix!(known_consts, :s, smallest_dim, 1)
        b.consts = known_consts
        idx = findfirst(==(smallest_dim), unknowns[:s])
        deleteat!(unknowns[:s], idx)
    end

    knowns = Channels{Vector{Field{T}}}(Tuple([
            V for V in fields(b.spectra[chan])
            if V in keys(b.consts[chan])
        ] for chan in chans
    ))

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(chans) - 1)
    nb_unknowns = [length(s) for s in b.spectra]
    nb_unknowns .-= length.([known_consts[chan] for chan in chans])

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
        LHS[1:nb_positions, col_idx] = .- facts[:t] .* b.block_values[:t][V]
        LHS[(nb_positions+1):end, col_idx] = zer
        col_idx += 1
    end
    for V in unknowns[:u]
        LHS[1:nb_positions, col_idx] = zer
        LHS[(nb_positions+1):end, col_idx] = .- facts[:u] .* b.block_values[:u][V]
        col_idx += 1
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = Vector{T}(undef, 2nb_positions)
    for i = 1:nb_positions
        RHS[i] = sum(b.block_values[:t][V][i] for V in knowns[:t]; init=zero(T))
        RHS[i] -= sum(b.block_values[:s][V][i] for V in knowns[:s]; init=zero(T))
        RHS[i+nb_positions] = sum(b.block_values[:u][V][i] for V in knowns[:u]; init=zero(T))
        RHS[i+nb_positions] -= sum(b.block_values[:s][V][i] for V in knowns[:s]; init=zero(T))
    end

    b.matrix = BootstrapMatrix{T}(unknowns, LHS, RHS)
end

function solve!(s::BootstrapSystem; precision_factor=1)
    # solve for two different sets of positions
    LHS, RHS = s.matrix.LHS, s.matrix.RHS
    if precision_factor > 1
        LHS = Matrix{Complex{BigFloat}}(LHS)
        RHS = Vector{Complex{BigFloat}}(RHS)
    end
    sol1, sol2 = setprecision(BigFloat, precision_factor * precision(BigFloat)) do
        (
            s.matrix.LHS[3:end, :] \ s.matrix.RHS[3:end],
            s.matrix.LHS[1:(end-2), :] \ s.matrix.RHS[1:(end-2)],
        )
    end

    # compare the two solutions
    errors = @. abs((sol1 - sol2) / sol2)

    # back to dictionary format
    mat = s.matrix
    chans = keys(s.spectra)
    nb_unknowns = vcat([0], [length(mat.unknowns[chan]) for chan in chans])
    for (i, chan) in enumerate(chans)
        for (j, V) in enumerate(mat.unknowns[chan])
            index = sum(nb_unknowns[k] for k = 1:(i-1); init=0) + j
            s.consts[chan][V] = sol1[index]
            s.consts.errors[chan][V] = errors[index]
        end
    end
end

function clear!(b::BootstrapSystem{T}) where {T}
    b.consts = StructureConstants(b.spectra)
    b.matrix = BootstrapMatrix{T}([chan for chan in keys(b.spectra)])
end
