function BootstrapMatrix{T}() where {T}
    BootstrapMatrix{T}(
        Channels{Vector{Tuple{Int,Field{T}}}}(
            Tuple(Vector{Field{T}}() for chan in (:s, :t, :u)),
        ),
        Matrix{T}(undef, 0, 0),
        Vector{T}(undef, 0),
    )
end

# function evaluate_blocks(bs::Vector{Block{T}}, xs) where {T}
#     nbs = length(bs)
#     results = Vector{Vector{T}}(undef, nbs)

#     Threads.@threads for i in 1:nbs
#         b = bs[i]
#         results[i] = [b(x) for x in xs]
#     end

#     return results
# end

function evaluate_blocks(bs::Vector{Block{T}}, xs) where {T}
    nthreads = Threads.nthreads()
    chunks = [bs[i:nthreads:end] for i = 1:nthreads]

    results_chunks = Vector{Vector{Vector{T}}}(undef, nthreads)

    Threads.@threads for tid = 1:nthreads
        local_chunk = chunks[tid]
        local_results = Vector{Vector{T}}(undef, length(local_chunk))
        for (j, b) in enumerate(local_chunk)
            local_results[j] = [b(x) for x in xs]
        end
        results_chunks[tid] = local_results
    end

    # Flatten result while preserving order
    return reduce(vcat, results_chunks)
end

function evaluate_blocks(S::Channels{U}, xs) where {T,U<:ChannelSpectrum{T}}
    return Channels{Vector{Vector{T}}}(
        Tuple(evaluate_blocks(S[chan].blocks, xs[chan]) for chan in keys(S)),
    )
end

function evaluate_blocks!(sys::BootstrapSystem)
    sys.blocks = evaluate_blocks(sys.spectra, sys.positions_cache)
    return sys.blocks
end

function new_random_point!(random_points, N; transfo = missing, square = true)
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
function random_points(N; transfo = missing, square = true)
    res = []
    for _ = 1:N
        new_random_point!(res, N, transfo = transfo, square = square)
    end
    return transfo !== missing ? transfo.(res) : res
end

function BootstrapSystem(
    S::Channels{U};
    knowns = nothing,
    extrapoints::Int = 6,
    pos = nothing,
) where {T,U<:ChannelSpectrum{T}}
    # if consts is empty, normalise the field with smallest indices
    chans = keys(S)
    co = S[1].corr
    if knowns === nothing
        knowns = StructureConstants{T}()
    end
    nb_unknowns = [length(s) for s in S]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])
    nb_lines = 2 * ((sum(nb_unknowns)+extrapoints) ÷ 2)
    if pos !== nothing
        nb_positions = length(pos)
    else
        nb_positions = nb_lines ÷ (length(S) - 1)
        pos = Vector{T}(random_points(nb_positions))
    end
    pos_chan = Channels{Vector{T}}(
        Tuple([channel_position(x, co, chan) for x in pos] for chan in keys(S)),
    )
    pos_cache = Channels{Vector{LRPositionCache{T}}}(
        Tuple([LRPositionCache(x, S[i].corr, chans[i]) for x in pos_chan[i]] for i = 1:3),
    )
    blocks = evaluate_blocks(S, pos_cache)
    matrix = BootstrapMatrix{T}()
    return BootstrapSystem{T,U}(pos, pos_cache, S, blocks, matrix, knowns)
end

function add!(b::BootstrapSystem, chan, V)
    new_block = add!(b.spectra[chan], V)
    if isnothing(new_block)
        error("trying to add already present block")
        return
    end
    # add a new position
    new_random_point!(b.positions, length(b.positions))
    new_pos = b.positions[end]
    for chan in (:s, :t, :u)
        new_pos_cache = LRPositionCache(new_pos, b.spectra[chan].corr, chan)
        push!(b.positions_cache[chan], new_pos_cache)
    end
    for chan in (:s, :t, :u)
        # evaluate previous blocks at the new position
        for v in b.blocks[chan]
            push!(v, new_block(new_pos))
        end
    end
    # evaluate the new block at all positions
    push!(b.blocks[chan], [new_block(x) for x in b.positions_cache[chan]])
end

function remove!(b::BootstrapSystem, chan, V)
    idx = remove!(b.spectra[chan], V)
    if isnothing(idx)
        error("the block you are trying to remove is not present")
        return
    end
    pop!(b.positions)
    for chan in (:s, :t, :u)
        pop!(b.positions_cache[chan])
    end
    # remove one position for all blocks
    for chan in (:s, :t, :u)
        for v in b.blocks[chan]
            pop!(v)
        end
    end
    # remove the block
    deleteat!(b.blocks[chan], idx);
end

function compute_linear_system!(b::BootstrapSystem{T}) where {T}
    chans = (:s, :t, :u)
    known_consts = b.consts

    unknowns = Channels{Vector{Tuple{Int,Field{T}}}}(
        Tuple(
            [
                (i, V) for (i, V) in enumerate(b.spectra[chan].fields) if
                !(V in keys(known_consts[chan]))
            ] for chan in chans
        ),
    )

    # if no consts are known, fix a normalisation: first str cst in the 1st-channel = 1
    if all([isempty(known_consts[chan]) for chan in chans])
        idx, V_norm = sort(unknowns[chans[1]], by = iV->real(total_dimension(iV[2])))[1]
        fix!(known_consts, chans[1], V_norm, 1)
        b.consts = known_consts
        deleteat!(unknowns[chans[1]], idx)
    end

    knowns = Channels{Vector{Tuple{Int,Field{T}}}}(
        Tuple(
            [
                (i, V) for (i, V) in enumerate(b.spectra[chan].fields) if
                V in keys(known_consts[chan])
            ] for chan in chans
        ),
    )

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(chans) - 1)
    nb_unknowns = [length(s) for s in b.spectra]
    nb_unknowns .-= length.([known_consts[chan] for chan in chans])

    # # solve Σ_unknowns (chan1 - chan2) = Σ_knowns (chan2 - chan1)
    # # and Σ_unknowns (chan1 - chan3) = Σ_knowns (chan2 - chan3)
    zer = zeros(T, nb_positions)

    # # Form the matrix
    # # [ G^s -G^t  0  ;
    # #   G^s  0   -G^u ]
    LHS = Matrix{T}(undef, nb_lines, sum(nb_unknowns))
    blocks = b.blocks[:s]
    _unknowns = unknowns[:s]
    for i = 1:length(_unknowns)
        LHS[1:nb_positions, i] = blocks[_unknowns[i][1]]
        LHS[(nb_positions+1):end, i] = blocks[_unknowns[i][1]]
    end
    blocks = b.blocks[:t]
    _unknowns = unknowns[:t]
    for i = 1:length(_unknowns)
        LHS[1:nb_positions, nb_unknowns[1]+i] = .- blocks[_unknowns[i][1]]
        LHS[(nb_positions+1):end, nb_unknowns[1]+i] = zer
    end
    blocks = b.blocks[:u]
    _unknowns = unknowns[:u]
    for i = 1:length(_unknowns)
        LHS[1:nb_positions, nb_unknowns[1]+nb_unknowns[2]+i] = zer
        LHS[(nb_positions+1):end, nb_unknowns[1]+nb_unknowns[2]+i] =
            .- blocks[_unknowns[i][1]]
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = Vector{T}(undef, 2nb_positions)
    for i = 1:nb_positions
        res = zero(T)
        for j = 1:length(knowns[:s])
            res -= b.blocks[:s][knowns[:s][j][1]][i]
        end
        for j = 1:length(knowns[:t])
            res += b.blocks[:t][knowns[:t][j][1]][i]
        end
        RHS[i] = res
        res = zero(T)
        for j = 1:length(knowns[:s])
            res -= b.blocks[:s][knowns[:s][j][1]][i]
        end
        for j = 1:length(knowns[:t])
            res += b.blocks[:t][knowns[:t][j][1]][i]
        end
        RHS[i+nb_positions] = res
    end
    # RHS = vcat(
    #     [sum(b.blocks[chans[i]][knowns[chans[i]][j][1]] for j in 1:length(knowns[chans[i]]); init=zer) .-
    #      sum(b.blocks[chans[1]][knowns[chans[i]][j][1]] for j in 1:length(knowns[chans[1]]); init=zer)
    #      for i in 2:length(chans)]...
    # )

    b.matrix = BootstrapMatrix{T}(unknowns, LHS, RHS)

    # b.matrix.LHS = LHS
    # b.matrix.RHS = RHS;
end

function solve!(s::BootstrapSystem; precision_factor = 1)
    @info "system size: $(size(s.matrix.LHS))"
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
        for (j, V) in mat.unknowns[chan]
            index = sum(nb_unknowns[k] for k = 1:(i-1); init = 0) + j
            s.consts[chan][V] = sol1[index]
            s.consts.errors[chan][V] = errors[index]
        end
    end
end

function clear!(b::BootstrapSystem{T}) where {T}
    b.consts = StructureConstants(b.spectra)
    b.matrix = BootstrapMatrix{T}([chan for chan in keys(b.spectra)]);
end
