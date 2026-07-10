"""
data required to build a system of bootstrap equations

# Fields

unknowns: list of fields whose structure constant is unknown, in a Channels struct
knowns: list of known (fields, structure_constant), in a Channels struct
moduli_count: number of moduli at which the conformal blocks are evaluated
block_values: values of the conformal blocks, in a Channels struct
factors: prefactors multiplying the block values. May depend on the moduli, the fields
rels: one of :st, :su, :tu, :stu. Imposes that the corresponding structure constants are
      equal, for instance rels = :st imposes that D^(s)_V = D^(t)_V.
"""
struct BootstrapSystem
    unknowns::Channels{Vector{CDorField}}
    knowns::Channels{Vector{Tuple{CDorField,Acb}}}
    moduli_count::Int
    block_values::Channels{Dict{CDorField,Vector{Acb}}}
    factors::Channels{Dict{CDorField,Vector{Acb}}}
    rels::Symbol
end

struct BootstrapMatrix
    fields
    A::Matrix{Acb}
    b::Vector{Acb}
end

function new_random_point(; square=true)
    xmin, xmax = square ? (0.1, 0.5) : (-0.4, 1.4)
    return z = xmin + (xmax - xmin) * rand() + (1 + 4 * rand()) * im / 10
end

function new_random_point!(random_points, N; square = true, cond = p -> true)
    while true
        z = new_random_point(square=square)
        sep = square ? 0.2 : 0.4
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
    Random.seed!(3456)
    res = []
    for _ = 1:N
        new_random_point!(res, N, square = square, cond = cond)
    end
    return transfo !== nothing ? transfo.(res) : res

    # Restore RNG state
    Random._GLOBAL_RNG[] = rng
end

choose_random_moduli(moduli_count, _::Correlation1) = 
    Vector{Acb}(random_points(moduli_count, transfo = x -> τfromx(x)))
choose_random_moduli(moduli_count, _::Correlation4) = 
    Vector{Acb}(random_points(moduli_count, transfo = nothing))

function BootstrapMatrix(s::BootstrapSystem)
    unknowns, knowns = s.unknowns, s.knowns
    moduli_count, block_values = s.moduli_count, s.block_values
    factors, rels = s.factors, s.rels
    if length(keys(block_values)) == 3
        if rels === :stu #= [(Gs-Gt);
                             (Gs-Gu)] =#
            cols_count = length(unknowns.s)
            offset_t = 0
            offset_u = 0
        elseif rels === :st #= [(Gs-Gt)  0 ;
                                 Gs     -Gu] =#

            offset_t = 0
            offset_u = length(unknowns.s)
            cols_count = length(unknowns.s) + length(unknowns.u)
        elseif rels === :su #= [Gs        -Gt;
                                (Gs-Gu)     0] =#
            offset_t = length(unknowns.s)
            offset_u = 0
            cols_count = length(unknowns.s) + length(unknowns.t)
        elseif rels === :tu #= [Gs     -Gt;
                                Gs     Gu)] =#
            offset_t = length(unknowns.s)
            offset_u = length(unknowns.s)
            cols_count = length(unknowns.s) + length(unknowns.t)
        elseif rels === nothing   #= [(Gs     -Gt        0
                                       Gs       0      -Gu] =#
            offset_t = length(unknowns.s)
            offset_u = length(unknowns.s) + length(unknowns.t)
            cols_count = length(unknowns.s) + length(unknowns.t) + length(unknowns.u)
        end

        A = zeros(Acb, 2*moduli_count, cols_count)
        b = zeros(Acb, 2*moduli_count)
        fields = Vector{Tuple{Symbol, CDorField}}(undef, cols_count)
        j = 0
        for chan in (:s, :t, :u)
            for V in unknowns[chan]
                j += 1
                fields[j] = (chan, V)
            end
        end
        for k in 1:moduli_count
            for (j, V) in enumerate(unknowns.s)
                A[k, j]                       += factors.s[V][k] * block_values.s[V][k]
                A[k+moduli_count, j]          += factors.s[V][k] * block_values.s[V][k]
            end
            for (j, V) in enumerate(unknowns.t)
                A[k, j+offset_t]              -= factors.t[V][k] * block_values.t[V][k]
            end
            for (j, V) in enumerate(unknowns.u)
                A[k+moduli_count, j+offset_u] -= factors.u[V][k] * block_values.u[V][k]
            end

            # right-hand side of the linear system
            for (V, D) in knowns.s
                b[k] -= factors.s[V][k] * D * block_values.s[V][k]
                b[k + moduli_count] -= factors.s[V][k] * D * block_values.s[V][k]
            end
            for (V, D) in knowns.t
                b[k] += factors.t[V][k] * D * block_values.t[V][k]
            end
            for (V, D) in knowns.u
                b[k + moduli_count] += factors.u[V][k] * D * block_values.u[V][k]
            end
        end
    else
        error("BootstrapSystem with 2 channels not implemented yet")
    end
    return BootstrapMatrix(unknowns, A, b)
end

function BootstrapSystem(co::Correlation, blocks::Channels;
                         fix=nothing, rels=nothing, moduli=nothing, extrapoints::Int=6,
                         prefactors=nothing)
    chans = collect(keys(blocks))
    chan_count = length(chans)

    unknowns_set = Channels(Dict(chan => Set(b.chan_field for b in blocks[chan])
                                 for chan in chans))
    knowns = Channels(Dict(chan => Tuple{CDorField, Acb}[] for chan in chans))

    # if nothing is fixed, fix the structure constant of the first found diagonal field to 1
    for chan in chans
        if fix !== nothing
            break
        end
        for V in unknowns_set[chan]
            if V.diagonal
                fix = [(chan, V, one(Acb))]
                break
            end
        end
    end
    # if no diagonal field is found, just fix the first structure constant we find in the s-channel to 1
    if fix === nothing
        V = first(unknowns_set.s)
        fix = [(:s, V, one(Acb))]
    end

    for (chan, V, D) in fix
        delete!(unknowns_set[chan], V)
        push!(knowns[chan], (V, D))
    end
    unknowns = Channels(Dict(chan => collect(unknowns_set[chan])
                             for chan in chans))
    
    unknowns_count = sum(length(unknowns[chan]) for chan in chans)
    lines_count = 2 * ((unknowns_count + extrapoints) ÷ 2)
    if moduli !== nothing
        moduli_count = length(moduli)
        if !(moduli_count > lines_count ÷ (chan_count - 1))
            error("You must give enough positions for the system to be overdetermined")
        end
    else
        moduli_count = chan_count == 1 ? lines_count ÷ 2 : lines_count ÷ (chan_count - 1)
    end

    moduli = choose_random_moduli(moduli_count, co)
     
    # compute the cache for all moduli, in each channel
    moduli_chan = Channels(Dict(chan => [channel_modulus(x, co, chan) for x in moduli]
                             for chan in chans))
    moduli_cache = Channels(Dict(chan => Vector(undef, length(moduli_chan[chan]))
                              for chan in chans)) # result container
    Threads.@threads for i = 1:chan_count
        chan = chans[i]
        moduli_cache[chan] = [PosCache(x, co, chan) for x in moduli_chan[chan]]
    end

    # evaluate the blocks
    thread_results = [Dict{CDorField,Vector{Acb}}() for _ = 1:max_thread_id()] # thread-local results
    block_values = Channels(Dict(chan => Dict{CDorField,Vector{Acb}}() for chan in chans))
    for chan in chans
        Threads.@threads for i = 1:length(blocks[chan])
            b = blocks[chan][i]
            thread_results[Threads.threadid()][b.chan_field] = [b(x) for x in moduli_cache[chan]]
        end
        # Merge all thread-local dictionaries into a single result
        for d in thread_results
            merge!(block_values[chan], d)
        end
    end

    # evaluate the prefactors
    if prefactors === nothing
        factors = Channels(Dict(chan => Dict(b.chan_field => ones(Acb, moduli_count)
                                             for b in blocks[chan]) for chan in chans))
    else
        factors = Channels(Dict(chan => Dict(b.chan_field => Vector(undef, moduli_count)
                                             for b in blocks[chan]) for chan in chans))
        for chan in chans
            for b in blocks[chan]
                for k in 1:moduli_count
                    factors[chan][V][k] = prefactors(chan, V, moduli[k])
                end
            end
        end
    end

    return BootstrapSystem(unknowns, knowns, moduli_count, block_values, factors, rels)
end

function solve!(s::BootstrapSystem)
    matrix = BootstrapMatrix(s)
    LHS = convert(Matrix{Complex{BigFloat}}, matrix.A)
    RHS = convert(Vector{Complex{BigFloat}}, matrix.b) # convert to BigFloat to use Julia's solver

    sol1 = LHS[3:end, :] \ RHS[3:end]
    sol2 = LHS[1:(end-2), :] \ RHS[1:(end-2)]

    # compare the two solutions
    errors = @. Float64(abs((sol1 - sol2) / sol2))

    res = Channels(Dict(chan => Dict() for chan in keys(s.block_values)))
    for (j, (chan, V)) in enumerate(matrix.fields)
        res[chan][V] = sol1[j]
    end
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

channel_modulus(x, _::Correlation4, chan) = crossratio(chan, x)
channel_modulus(x, _::Correlation1, chan) = modular_param(chan, x)
