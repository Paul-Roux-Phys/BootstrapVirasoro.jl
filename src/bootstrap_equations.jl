struct BootstrapMatrix{T}
    unknowns::Channels{Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
end

"""
TODO
"""
mutable struct StructureConstants{T}
    constants::Channels{Dict{Field{T}, T}}
    errors::Channels{Dict{Field{T}, T}}
end

mutable struct BootstrapSystem{T, U<:ChannelSpectrum{T}}
    positions::Vector{T}
    positions_cache::Channels{Vector{LRPositionCache{T}}}        # positions at which eqs are evaluated
    spectra::Channels{U}    # channel spectra
    blocks::Channels{Vector{Vector{T}}} # all blocks evaluated
    matrix::BootstrapMatrix{T}  # matrix of equations
    consts::StructureConstants{T}
end

function BootstrapMatrix{T}() where {T}
    BootstrapMatrix{T}(
        Channels{Vector{Field{T}}}(Tuple(
            Vector{Field{T}}()
            for chan in (:s, :t, :u)
        )),
        Matrix{T}(undef, 0, 0),
        Vector{T}(undef, 0)
    )
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
    xmin, xmax, sep = square ? (.1, .5, .2) : (-.4, 1.4, .4)
    res = []
    while length(res) < N
        z = xmin + (xmax - xmin)*rand() + (1+4*rand())*im/10 
        if minimum(append!(abs.(z .- res), 1)) > sep/sqrt(N)
            append!(res, z)
        end
    end
    return transfo !== missing ? transfo.(res) : res
end

function evaluate_blocks(bs::Vector{Block{T}}, xs) where {T}
    nbs = length(bs)
    results = Vector{Vector{T}}(undef, nbs)

    Threads.@threads for i in 1:nbs
        b = bs[i]
        results[i] = [b(x) for x in xs]
    end

    return results
end

function evaluate_blocks(S::Channels{U}, xs) where {T, U<:ChannelSpectrum{T}}
    return Channels{Vector{Vector{T}}}(Tuple(
        evaluate_blocks(S[chan].blocks, xs[chan])
        for chan in keys(S)
    ))
end

function BootstrapSystem(S::Channels{U}; knowns=nothing, extrapoints::Int=6) where {T, U<:ChannelSpectrum{T}}
    # if consts is empty, normalise the field with smallest indices
    chans = keys(S)
    co = S[1].corr
    if knowns === nothing
        knowns = StructureConstants{T}()
    end
    nb_unknowns = [length(s) for s in S]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])
    nb_lines = 2 * ((sum(nb_unknowns)+extrapoints) ÷ 2)
    nb_positions = nb_lines ÷ (length(S) - 1)
    pos = Vector{T}(random_points(nb_positions))
    pos_chan = Channels{Vector{T}}(Tuple(
        [channel_position(x, co, chan) for x in pos]
        for chan in keys(S)
    ))
    pos_cache = Channels{Vector{LRPositionCache{T}}}(Tuple(
        [LRPositionCache(x, S[i].corr, chans[i]) for x in pos_chan[i]]
        for i in 1:3
    ))
    blocks = evaluate_blocks(S, pos_cache)
    matrix = BootstrapMatrix{T}()
    return BootstrapSystem{T, U}(
        pos, pos_cache, S, blocks, matrix, knowns
    )
end

function compute_linear_system!(b::BootstrapSystem{T}) where {T}
    chans = (:s, :t, :u)
    knowns = b.consts

    # if no consts are known, fix a normalisation: first str cst in the 1st-channel = 1
    unknowns = Channels{Vector{Field{T}}}(Tuple(
        b.spectra[chan].fields for chan in chans
    ))
    if all([isempty(b.consts[chan]) for chan in chans])
        idx, V_norm = sort(
            collect(enumerate(unknowns[chans[1]])),
            by=V->real(total_dimension(V[2]))
        )[1]
        fix!(knowns, chans[1], V_norm, 1)
        b.consts = knowns
        deleteat!(unknowns[chans[1]], idx)
    end

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(chans) - 1)
    nb_unknowns = [length(s) for s in b.spectra]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])

    # # solve Σ_unknowns (chan1 - chan2) = Σ_knowns (chan2 - chan1)
    # # and Σ_unknowns (chan1 - chan3) = Σ_knowns (chan2 - chan3)
    zer = zeros(T, nb_positions)

    # # Form the matrix
    # # [ G^s -G^t  0  ;
    # #   G^s  0   -G^u ]
    LHS = Matrix{T}(undef, nb_lines, sum(nb_unknowns))
    for i in 1:length(unknowns[chans[1]])
        LHS[1:nb_positions, i] = b.blocks[chans[1]][i]
        LHS[nb_positions+1:end, i] = b.blocks[chans[1]][i]
    end
    for i in 1:length(unknowns[chans[2]])
        LHS[1:nb_positions, nb_unknowns[1]+i] = .- b.blocks[chans[2]][i]
        LHS[nb_positions+1:end, nb_unknowns[1]+i] = zer
    end
    for i in 1:length(unknowns[chans[3]])
        LHS[1:nb_positions, nb_unknowns[1]+nb_unknowns[2]+i] = zer
        LHS[nb_positions+1:end, nb_unknowns[1]+nb_unknowns[2]+i] = .- b.blocks[chans[3]][i]
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = vcat(
        [sum(b.blocks[chans[i]][j] for j in 1:length(knowns[chans[i]]); init=zer) .-
         sum(b.blocks[chans[1]][j] for j in 1:length(knowns[chans[1]]); init=zer)
         for i in 2:length(chans)]...
    )

    b.matrix = BootstrapMatrix(unknowns, LHS, RHS)

    # b.matrix.LHS = LHS
    # b.matrix.RHS = RHS;
end

function solve!(s::BootstrapSystem; precision_factor=1)
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
            s.matrix.LHS[1:end-2, :] \ s.matrix.RHS[1:end-2]
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
            index = sum(nb_unknowns[k] for k in 1:i-1; init=0) + j
            s.consts[chan][V] = sol1[index]
            s.consts.errors[chan][V] = errors[index]
        end
    end
end

function clear!(b::BootstrapSystem{T}) where {T}
    b.consts = StructureConstants(b.spectra)
    b.matrix = BootstrapMatrix{T}([chan for chan in keys(b.spectra)]);
end
