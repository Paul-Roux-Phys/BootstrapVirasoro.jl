# fit a polynomial in one variable X
# degs: the degree
# data: the numerical data to fit: pairs of (x, value)
function fit1(deg, data)
    mat = [x^n for (x, _) in data, n = 0:deg]
    vals = [v for (_, v) in data]
    mat \ vals
end

data1 = [(0, 1), (1, 2), (2, 3), (3, 4)]
fit1(1, data1)

function fit2(deg1, deg2, data)
    @assert length(data) > (deg1 + 1) * (deg2 + 1) "
You must provide more data points
"
    mat = Matrix{float(typeof(data).parameters[1].parameters[2])}(
        undef,
        length(data),
        (deg1 + 1) * (deg2 + 1),
    )
    monoms = [(m, n) for m = 0:deg1 for n = 0:deg2]
    for i in eachindex(data)
        x, y = data[i][1]
        for j in eachindex(monoms)
            mat[i, j] = x^(monoms[j][1]) * y^(monoms[j][2])
        end
    end
    vals = [v for (_, v) in data]
    mat \ vals
end

data2 = [((x, 0), x + 1) for x = 1:10]
fit2(2, 2, data2)

function fit(degs, data)
    nb_monoms = prod(d + 1 for d in degs)
    @assert length(data) > nb_monoms "
You must provide more data points
"
    T = float(typeof(data).parameters[1].parameters[2])
    mat = Matrix{T}(undef, length(data), nb_monoms)
    ranges = (0:d for d in degs)
    monoms = [t for t in Iterators.product(ranges...)]
    for (i, d) in enumerate(data)
        xvec = d[1]
        for (j, m) in enumerate(monoms)
            mat[i, j] = prod(x^m[k] for (k, x) in enumerate(xvec))
        end
    end
    vals = [v for (_, v) in data]
    res = mat \ vals
    return [(monoms[i], res[i]) for i in eachindex(monoms)]
end

data = [
    ((x, y), 2x^2 + 3x*y +y) for x = 1:0.4:4 for y in 0:0.3:0.8
]
fit([2, 1], data)
