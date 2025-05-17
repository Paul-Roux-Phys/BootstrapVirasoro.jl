struct BootstrapMatrix{T}
    channels::Vector{Symbol}
    fields::Dict{Symbol,Vector{Field{T}}}
    LHS::Matrix{T}
    RHS::Vector{T}
end

function BootstrapMatrix{T}(chans) where {T}
    BootstrapMatrix{T}(
        chans,
        Dict(chan => Vector{Field{T}}[] for chan in chans),
        Matrix{T}[],
        Vector{T}[]
    )
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
    blocks::Dict{Symbol, Vector{Block{T}}} # all blocks evaluated
    matrix::BootstrapMatrix{T}  # matrix of equations, possibly with excluded fields
    consts::StructureConstants{T}
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

function BootstrapMatrix(T, chans, knowns, unknowns, nb_unknowns, nb_positions, nb_lines, blocks)
    return BootstrapMatrix{T}(chans, unknowns, matrix, RHS)
end

function BootstrapSystem(S::Dict{Symbol, U}, consts; extrapoints::Int=6) where {T, U<:ChannelSpectrum{T}}
    # if consts is empty, normalise the field with smallest indices
    chans = [chan for chan in keys(S)]
    nb_unknowns = [length(s) for (chan, s) in S]
    nb_unknowns .-= length.([consts[chan] for chan in chans])
    nb_lines = 2 * ((sum(nb_unknowns)+extrapoints) ÷ 2)
    nb_positions = nb_lines ÷ (length(S) - 1)
    positions = Vector{T}(random_points(nb_positions))
    blocks = evaluate_blocks(S, positions)
    # knowns = [consts[chan] for chan in chans]

    # matrix = BootstrapMatrix(
    #     T, chans, knowns, unknowns, nb_unknowns, nb_positions, nb_lines, blocks
    # )
    matrix = BootstrapMatrix{T}(chans)
    return BootstrapSystem{T, U}(positions, S, blocks, matrix, consts)
end

function compute_matrix!(b::BootstrapSystem{T}; excludes=[], knowns=nothing) where {T}
    chans = b.matrix.channels
    unknowns = Dict(
        chan => [
            V for V in keys(b.blocks[chan]) if !(V in excludes) && !(V in b.consts.fields)
        ]
        for chan in chans
    )

    if knowns === nothing
        # if no consts are known, fix a normalisation: first str cst in the s-channel = 1
        V_norm = sort(
            b.spectra[:s].fields,
            by=V->real(total_dimension(V))
        )[1]
        knowns = StructureConstants(b.spectra)
        knowns[:s][V_norm] = 1
    end
    
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
        [sum(b.blocks[chans[i]][V] for V in keys(knowns[i]); init=zer) .-
         sum(b.blocks[chans[1]][V] for V in keys(knowns[1]); init=zer)
         for i in 2:length(chans)]...
    )
    b.matrix.LHS = LHS
    b.matrix.RHS = RHS
end
