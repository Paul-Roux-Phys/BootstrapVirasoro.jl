"""
TODO
"""
mutable struct StructureConstants{T}
    constants::Channels{Dict{Field{T}, T}}
    errors::Channels{Dict{Field{T}, T}}
end

function StructureConstants{T}() where {T}
    constants = Channels{Dict{Field{T}, T}}(Tuple(
        Dict() for chan in (:s, :t, :u)
    ))
    errors = deepcopy(constants)
    return StructureConstants{T}(constants, errors)
end

function Base.getproperty(c::StructureConstants, s::Symbol)
    s === :fields && begin
        consts = getfield(c, :constants)
        return vcat([[V for V in keys(consts[chan])] for chan in keys(consts)]...)
    end
    getfield(c, s)
end

Base.getindex(c::StructureConstants, s::Symbol) = c.constants[s]

function fix!(cnst, chan, field, value; error = 0)
    cnst[chan][field] = value
    cnst.errors[chan][field] = error
end

struct BootstrapMatrix{T}
    unknowns::Channels{Vector{Tuple{Int, Field{T}}}}
    LHS::Matrix{T}
    RHS::Vector{T}
end

function BootstrapMatrix{T}() where {T}
    BootstrapMatrix{T}(
        Channels{Vector{Tuple{Int, Field{T}}}}(Tuple(
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

mutable struct BootstrapSystem{T, U<:ChannelSpectrum{T}}
    positions::Vector{T}
    positions_cache::Channels{Vector{LRPositionCache{T}}}        # positions at which eqs are evaluated
    spectra::Channels{U}    # channel spectra
    blocks::Channels{Vector{Vector{T}}} # all blocks evaluated
    matrix::BootstrapMatrix{T}  # matrix of equations
    consts::StructureConstants{T}
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
    known_consts = b.consts

    unknowns = Channels{Vector{Tuple{Int, Field{T}}}}(Tuple(
        [(i, V) for (i, V) in enumerate(b.spectra[chan].fields)
             if !(V in keys(known_consts[chan]))]
        for chan in chans
    ))

    # if no consts are known, fix a normalisation: first str cst in the 1st-channel = 1
    if all([isempty(known_consts[chan]) for chan in chans])
        idx, V_norm = sort(
            unknowns[chans[1]],
            by=iV->real(total_dimension(iV[2]))
        )[1]
        fix!(known_consts, chans[1], V_norm, 1)
        b.consts = known_consts
        deleteat!(unknowns[chans[1]], idx)
    end

    knowns = Channels{Vector{Tuple{Int, Field{T}}}}(Tuple(
        [(i, V) for (i, V) in enumerate(b.spectra[chan].fields)
             if V in keys(known_consts[chan])]
        for chan in chans
    ))

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
    for i in 1:length(_unknowns)
        LHS[1:nb_positions, i] = blocks[_unknowns[i][1]]
        LHS[nb_positions+1:end, i] = blocks[_unknowns[i][1]]
    end
    blocks = b.blocks[:t]
    _unknowns = unknowns[:t]
    for i in 1:length(_unknowns)
        LHS[1:nb_positions, nb_unknowns[1]+i] = .- blocks[_unknowns[i][1]]
        LHS[nb_positions+1:end, nb_unknowns[1]+i] = zer
    end
    blocks = b.blocks[:u]
    _unknowns = unknowns[:u]
    for i in 1:length(_unknowns)
        LHS[1:nb_positions, nb_unknowns[1]+nb_unknowns[2]+i] = zer
        LHS[nb_positions+1:end, nb_unknowns[1]+nb_unknowns[2]+i] = .- blocks[_unknowns[i][1]]
    end

    # Form the right hand side: [  Σ_knowns (chan2 - chan1),  Σ_knowns (chan3 - chan1) ]
    RHS = Vector{T}(undef, 2nb_positions)
    for i in 1:nb_positions
        res = zero(T)
        for j in 1:length(knowns[:s])
            res -= b.blocks[:s][knowns[:s][j][1]][i]
        end
        for j in 1:length(knowns[:t])
            res += b.blocks[:t][knowns[:t][j][1]][i]
        end
        RHS[i] = res
        res = zero(T)
        for j in 1:length(knowns[:s])
            res -= b.blocks[:s][knowns[:s][j][1]][i]
        end
        for j in 1:length(knowns[:t])
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
        for (j, V) in mat.unknowns[chan]
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

function format_complex(z::Complex{<:Real})
    real_str = @sprintf("%.5e", real(z))
    imag_str = @sprintf("%.5e", abs(imag(z)))
    sign = imag(z) < 0 ? "-" : "+"
    buf = real(z) > 0 ? " " : ""
    return "$buf$real_str $sign $(imag_str)im"
end

function Base.show(io::IO, c::StructureConstants{T}) where {T}
    chans = keys(c.constants)
    nondiags = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if !isdiagonal(V)],
            by=V -> (V.r, V.s)
        )
        for chan in chans
    )
    diags = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if isdiagonal(V) && !isdegenerate(V)],
            by=V -> real(total_dimension(V))
        )
        for chan in chans
    )
    degs = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if isdegenerate(V)],
            by=V -> real(total_dimension(V))
        )
        for chan in chans
    )

    for chan in sort(collect(chans), by=string)
        # Collect all the labels for this channel
        all_Vs = vcat(degs[chan], diags[chan], nondiags[chan])

        # Find max width of label for alignment
        max_label_width = if isempty(all_Vs)
            0
        else
            maximum(length(string(V)) for V in all_Vs)
        end

        str_cst_col_width = if isempty(all_Vs)
            0
        else
            length(format_complex(zero(T)))
        end

        # Print channel header
        println(io, "Channel $(chan)")
        println(io, repeat('=', 9 + length(string(chan))))
        
        # Print column headers
        label1 = rpad("Fields", max_label_width)
        label2 = rpad("Structure constants", str_cst_col_width)
        label3 = rpad("Relative errors", str_cst_col_width)
        println(io, "$label1 | $label2  | $label3")
        println(io, repeat("-", max_label_width+2*str_cst_col_width))

        # Print each row with alignment
        for V in all_Vs
            label = rpad(string(V), max_label_width)
            value = format_complex(c[chan][V])
            error = format_complex(c.errors[chan][V])
            println(io, "$label | $value | $error")
        end
    end
end

Base.length(c::StructureConstants) = sum(length(c.constants[chan]) for chan in keys(c.constants))
