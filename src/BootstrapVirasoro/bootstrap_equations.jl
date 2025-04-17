export crossing_equations

const Spectra{T} = Tuple{Vararg{Spectrum{T}}}

"""
Generates N points in the square (.1, .5) + (.1, .5)i, while respecting 
a minimum distance .2/sqrt(N) between points. Or in the rectangle 
(-.4, 1.4) + (.1, .5)i with a distance .4/sqrt(N).
N = the number of points.
function = a function that we may apply on the points
square = whether to use a square (otherwise, a rectangle)
"""
function draw_points(N; transfo=missing, square=true)
    xmin, xmax, sep = square ? (.1, .5, .2) : (-.4, 1.4, .4)
    res = []
    while length(res) < N
        z = xmin + (xmax - xmin)*rand() + (1+4*rand())*im/10 
        if minimum(append!(abs.(z .- res), 1)) > sep/sqrt(N)
            append!(res, z)
        end
    end
    transfo !== missing ? transfo.(res) : res
end

function crossing_equations(co, S::Spectra{T}, x::Complex) where {T}
    if length(S) == 2
        return hcat(
            [evaluate(b, x) for b in S[1].blocks]',
            [-evaluate(b, x) for b in S[2].blocks]'
        )
    else # length(S) == 3
        s_chan = [evaluate(b, x) for b in S[1].blocks]'
        return vcat(
            hcat(
                s_chan,
                [-evaluate(b, x) for b in S[2].blocks]',
                zeros(size(S[3].blocks))'
            ),
            hcat(
                s_chan,
                zeros(size(S[2].blocks))',
                [-evaluate(b, x) for b in S[3].blocks]'
            )
        )
    end
end

function crossing_equations(co, S::Spectra{T}; extrapoints::Int=5) where {T}
    nb_unknowns = sum(length(s) for s in S)
    nb_positions = (nb_unknowns+extrapoints) ÷ (length(S) - 1)
    positions = draw_points(nb_positions)
    vcat([crossing_equations(co, S, x) for x in positions]...)
end

function create_matrix(
    co::Correlation,
    S::Vector{Spectrum};
    extrapoints=5,
    knowns = 1
)
    # count unknowns
    N_unknowns = sum(length(s) for (c, s) in S) - knowns
    N_channels = length(s)
    N_points = (N_unknowns + extrapoints) ÷ (N_channels - 1)

    # draw positions
    positions = 2
end
