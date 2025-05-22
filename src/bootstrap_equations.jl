mutable struct BootstrapMatrix{T}
    channels::Vector{Symbol}
    unknowns::Dict{Symbol,Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
end

"""
TODO
"""
mutable struct StructureConstants{T}
    constants::Dict{Symbol,Dict{Field{T},T}}
    errors::Dict{Symbol,Dict{Field{T},T}}
end

mutable struct BootstrapSystem{T, U<:ChannelSpectrum{T}}
    positions::Vector{T}        # positions at which eqs are evaluated
    spectra::Dict{Symbol, U}    # channel spectra
    blocks::Dict{Symbol, Dict{Field{T}, Vector{T}}} # all blocks evaluated
    matrix::BootstrapMatrix{T}  # matrix of equations
    consts::StructureConstants{T}
end

function BootstrapMatrix{T}(chans) where {T}
    BootstrapMatrix{T}(
        chans,
        Dict(chan => Vector{Field{T}}[] for chan in chans),
        Matrix{T}(undef, 0, 0),
        Vector{T}(undef, 0)
    )
end

Base.getindex(c::StructureConstants, s::Symbol) = c.constants[s]

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

function evaluate_blocks(S, xs)
    blocks = Dict(
        chan => Dict(b.channel_field => b(xs) for b in s.blocks)
        for (chan, s) in S
    )
    return blocks
end

function BootstrapSystem(S::Dict{Symbol, U}; knowns=nothing, extrapoints::Int=6) where {T, U<:ChannelSpectrum{T}}
    # if consts is empty, normalise the field with smallest indices
    chans = [chan for chan in keys(S)]
    if knowns === nothing
        knowns = StructureConstants{T}(chans)
    end
    nb_unknowns = [length(s) for (chan, s) in S]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])
    nb_lines = 2 * ((sum(nb_unknowns)+extrapoints) ÷ 2)
    nb_positions = nb_lines ÷ (length(S) - 1)
    positions = Vector{T}(random_points(nb_positions))
    blocks = evaluate_blocks(S, positions)
    matrix = BootstrapMatrix{T}(chans)
    return BootstrapSystem{T, U}(positions, S, blocks, matrix, knowns)
end

function compute_linear_system!(b::BootstrapSystem{T}) where {T}
    chans = b.matrix.channels
    knowns = b.consts

    # if no consts are known, fix a normalisation: first str cst in the 1st-channel = 1
    unknowns = Dict(
        chan => [
            V for V in keys(b.blocks[chan]) if !(V in keys(knowns[chan]))
        ]
        for chan in chans
    )
    if all([isempty(b.consts[chan]) for chan in chans])
        idx, V_norm = sort(
            collect(enumerate(unknowns[chans[1]])),
            by=V->real(total_dimension(V[2]))
        )[1]
        fix!(knowns, chans[1], V_norm, 1)
        b.consts = knowns
        deleteat!(unknowns[chans[1]], idx)
    end
    b.matrix.unknowns = unknowns

    nb_positions = length(b.positions)
    nb_lines = nb_positions * (length(chans) - 1)
    nb_unknowns = [length(s) for (chan, s) in b.spectra]
    nb_unknowns .-= length.([knowns[chan] for chan in chans])

    # solve Σ_unknowns (chan1 - chan2) = Σ_knowns (chan2 - chan1)
    # and Σ_unknowns (chan1 - chan3) = Σ_knowns (chan2 - chan3)
    zer = zeros(T, nb_positions)

    # Form the matrix
    # [ G^s -G^t  0  ;
    #   G^s  0   -G^u ]
    LHS = Matrix{T}(undef, nb_lines, sum(nb_unknowns))
    for (i, V) in enumerate(unknowns[chans[1]])
        LHS[1:nb_positions, i] = b.blocks[chans[1]][V]
        LHS[nb_positions+1:end, i] = b.blocks[chans[1]][V]
    end
    for (i, V) in enumerate(unknowns[chans[2]])
        LHS[1:nb_positions, nb_unknowns[1]+i] = b.blocks[chans[2]][V]
        LHS[nb_positions+1:end, nb_unknowns[1]+i] = zer
    end
    for (i, V) in enumerate(unknowns[chans[3]])
        LHS[1:nb_positions, nb_unknowns[1]+nb_unknowns[2]+i] = zer
        LHS[nb_positions+1:end, nb_unknowns[1]+nb_unknowns[2]+i] = b.blocks[chans[3]][V]
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = vcat(
        [sum(b.blocks[chans[i]][V] for V in keys(knowns[chans[i]]); init=zer) .-
         sum(b.blocks[chans[1]][V] for V in keys(knowns[chans[1]]); init=zer)
         for i in 2:length(chans)]...
    )

    b.matrix.LHS = LHS
    b.matrix.RHS = RHS;
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
    nb_unknowns = vcat([0], [length(mat.unknowns[chan]) for chan in mat.channels])
    for (i, chan) in enumerate(mat.channels)
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
