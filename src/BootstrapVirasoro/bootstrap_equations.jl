const Spectra = Dict{Symbol, BulkSpectrum};

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
        z = xmin + (xmax - xmin)*rand() + (1+4*rand())*1j/10 
        if min([abs(z - point) for point in points] + [1]) > sep/sqrt(N)
            append!(res, transfo !== missing ? transfo(z) : z)
        end
    end

    res
end

function create_matrix(
    c::Correlation,
    S::Spectra;
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