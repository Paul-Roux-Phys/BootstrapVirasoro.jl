export BootstrapSystem
using LinearAlgebra: I, norm
using GenericLinearAlgebra: qr
export StructureConstants

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
        Matrix{T}[]
        Vector{T}[]
    )
end

"""
Structure constants associated to a correlation function
Obtained by solving bootstrap equations.
"""
mutable struct StructureConstants{T}
    constants::Dict{Symbol,Dict{Field{T},T}}
    errors::Dict{Symbol,Dict{Field{T},T}}
end

mutable struct BootstrapSystem{T, U<:ChannelSpectrum{T}}
    positions::Vector{T}        # positions at which eqs are evaluated
    spectra::Dict{Symbol, U}    # channel spectra
    matrix::BootstrapMatrix{T}
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
    # solve Σ_unknowns (chan1 - chan2) = Σ_knowns (chan2 - chan1)
    # and Σ_unknowns (chan1 - chan3) = Σ_knowns (chan2 - chan3)
    zer = zeros(T, nb_positions)

    # Form the matrix
    # [ G^s -G^t  0  ;
    #   G^s  0   -G^u ]
    matrix = Matrix{T}(undef, nb_lines, sum(nb_unknowns))
    for (i, V) in enumerate(unknowns[chans[1]])
        matrix[1:nb_positions, i] = blocks[chans[1]][V]
        matrix[nb_positions+1:end, i] = blocks[chans[1]][V]
    end
    for (i, V) in enumerate(unknowns[chans[2]])
        matrix[1:nb_positions, nb_unknowns[1]+i] = blocks[chans[2]][V]
        matrix[nb_positions+1:end, nb_unknowns[1]+i] = zer
    end
    for (i, V) in enumerate(unknowns[chans[3]])
        matrix[1:nb_positions, nb_unknowns[1]+nb_unknowns[2]+i] = zer
        matrix[nb_positions+1:end, nb_unknowns[1]+nb_unknowns[2]+i] = blocks[chans[3]][V]
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = vcat(
        [sum(blocks[chans[i]][V] for V in keys(knowns[i]); init=zer) .-
         sum(blocks[chans[1]][V] for V in keys(knowns[1]); init=zer)
         for i in 2:length(chans)]...
    )
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
    unknowns = Dict(
        chan => [field for field in keys(blocks[chan]) if !(field in keys(consts[chan]))]
        for chan in chans
    )
    knowns = [consts[chan] for chan in chans]

    matrix = BootstrapMatrix(
        T, chans, knowns, unknowns, nb_unknowns, nb_positions, nb_lines, blocks
    )
    return BootstrapSystem{T, U}(positions, S, matrix, consts)
end

function BootstrapSystem(S, consts; extrapoints::Int=6, exclude=nothing)
    
end
