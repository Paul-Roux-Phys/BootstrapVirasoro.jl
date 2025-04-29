export bootstrap_equations

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

function bootstrap_equations(S, x)
    if use_distributed()
        return Dict(pmap(
            (chan, s) -> (chan => [b(x) for b in s.blocks]),
            collect(S)
        ))
    end
    blocks = Dict(
        chan => [b(x) for b in s.blocks]
        for (chan, s) in S
    )
    # if length(S) == 2
    #     return hcat(res[1]', .-res[2]')
    # elseif length(S) == 3
    #     return vcat(
    #         hcat(res[1]', .-res[2]', zeros(size(S[3].blocks))'),
    #         hcat(res[1]', zeros(size(S[2].blocks))', .-res[3]')
    #     ) 
    # else error("You can only input 2 or 3 channels")
    # end
end

function bootstrap_equations(co, S::Spectra{T}; extrapoints::Int=5) where {T}
    nb_unknowns = sum(length(s) for s in S)
    nb_positions = (nb_unknowns+extrapoints) ÷ (length(S) - 1)
    positions = draw_points(nb_positions)
    bootstrap_equations(co, S, positions)
end
