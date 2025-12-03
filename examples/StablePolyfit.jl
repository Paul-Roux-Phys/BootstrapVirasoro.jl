using LinearAlgebra, BootstrapVirasoro, Printf
using LinearAlgebra, Printf

# ---------------------------------------------------------
# Chebyshev polynomials (first kind)
# ---------------------------------------------------------
function chebyshevT(k::Int, x)
    k == 0 && return 1
    k == 1 && return x
    T0, T1 = 1, x
    for _ in 2:k
        T0, T1 = T1, 2*x*T1 - T0
    end
    return T1
end

# ---------------------------------------------------------
# Generate exponent tuples for total degree ≤ D
# ---------------------------------------------------------
function total_degree_exponents(nvars::Int, D::Int)
    exps = Tuple[]
    function rec(prefix, var, sumdeg)
        var > nvars && return push!(exps, Tuple(prefix))
        for d in 0:(D - sumdeg)
            rec([prefix...; d], var + 1, sumdeg + d)
        end
    end
    rec(Int[], 1, 0)
    return exps
end

# ---------------------------------------------------------
# Main type
# ---------------------------------------------------------
mutable struct Polyfit
    varnames::Vector{Symbol}
    maxdegree::Int
    exponents::Vector{Tuple{Vararg{Int}}}
    scale_min::Vector{Float64}
    scale_max::Vector{Float64}
    scale_center::Vector{Float64}  # automatically center
    coeffs::Vector{Float64}
end

# ---------------------------------------------------------
# Constructor
# ---------------------------------------------------------
function Polyfit(varnames::Vector{Symbol}, maxdegree::Int)
    n = length(varnames)
    exps = total_degree_exponents(n, maxdegree)
    return Polyfit(varnames, maxdegree, exps,
                   Float64[], Float64[], Float64[], Float64[])
end

# ---------------------------------------------------------
# Scaling utilities
# ---------------------------------------------------------
scale_to_cheb(x, c, s) = 2*(x-c)/s  # maps [c-s/2, c+s/2] -> [-1,1]

function scale_point!(pf::Polyfit, x::Vector)
    [scale_to_cheb(x[i], pf.scale_center[i], pf.scale_max[i] - pf.scale_min[i])
        for i in eachindex(x)]
end

# ---------------------------------------------------------
# Build design matrix row
# ---------------------------------------------------------
function basis_row(pf::Polyfit, x_scaled)
    n = length(x_scaled)
    maxdeg = pf.maxdegree
    # precompute Chebyshev values
    Tvals = [ [chebyshevT(k, x_scaled[i]) for k in 0:maxdeg] for i in 1:n ]
    # tensor product
    row = similar(pf.exponents, Float64)
    for (j, e) in enumerate(pf.exponents)
        prod_val = 1.0
        @inbounds for i in 1:n
            prod_val *= Tvals[i][e[i]+1]
        end
        row[j] = prod_val
    end
    return row
end

# ---------------------------------------------------------
# Fit routine
# ---------------------------------------------------------
function fit!(pf::Polyfit, data)
    nvars = length(pf.varnames)

    # extract X and Y
    X = [d[1] for d in data]
    Y = [d[2] for d in data]

    # compute scaling and centering
    pf.scale_min = [minimum(x[i] for x in X) for i in 1:nvars]
    pf.scale_max = [maximum(x[i] for x in X) for i in 1:nvars]
    pf.scale_center = [(pf.scale_min[i] + pf.scale_max[i])/2 for i in 1:nvars]

    # build design matrix
    M = Matrix{Float64}(undef, length(X), length(pf.exponents))
    for i in 1:length(X)
        xscaled = scale_point!(pf, X[i])
        M[i, :] = basis_row(pf, xscaled)
    end

    # QR-based LSQ
    pf.coeffs = qr(M) \ Y
    return pf
end

# ---------------------------------------------------------
# LaTeX printing
# ---------------------------------------------------------
import Base: show

function show(io::IO, ::MIME"text/latex", pf::Polyfit)
    terms = String[]
    for (c, e) in zip(pf.coeffs, pf.exponents)
        term = string(
            @sprintf("%.5f", c), " ",
            join(["T_{$(e[i])}($(pf.varnames[i]))" for i in 1:length(e)], " ")
        )
        push!(terms, term)
    end
    print(io, join(terms, " + "))
end

# ---------------------------------------------------------
# Convert a univariate Chebyshev polynomial to monomial coeffs
# Returns a Vector{Float64} where index = power + 1
# e.g. [c0, c1, c2] represents c0 + c1*x + c2*x^2
# ---------------------------------------------------------
function cheb_to_monomial_coeffs(n)
    if n == 0
        return [1.0]
    elseif n == 1
        return [0.0, 1.0]
    else
        Tnm2 = [1.0]       # T0
        Tnm1 = [0.0, 1.0]  # T1
        for k in 2:n
            # 2*x*T_{n-1} - T_{n-2}
            # multiply Tnm1 by 2x: shift by 1 and multiply by 2
            shifted = [0.0; 2.0 .* Tnm1]
            # subtract Tnm2, pad Tnm2 with zeros
            len = max(length(shifted), length(Tnm2))
            coeffs = [ (i <= length(shifted) ? shifted[i] : 0.0) -
                       (i <= length(Tnm2) ? Tnm2[i] : 0.0)
                       for i in 1:len ]
            Tnm2, Tnm1 = Tnm1, coeffs
        end
        return Tnm1
    end
end
"""
    polyfit_to_monomial(pf::StablePolyfit.Polyfit)

Convert a multivariate Chebyshev polynomial fit to the **monomial basis in original x coordinates**.
Handles arbitrary intervals and automatically accounts for centering and scaling.
Returns a Dict mapping exponent tuples `(i,j,...) => coefficient`.
"""
function polyfit_to_monomial(pf::Polyfit)
    monomial_coeffs = Dict{Tuple{Vararg{Int}}, Float64}()

    nvars = length(pf.varnames)
    s = [pf.scale_max[i] - pf.scale_min[i] for i in 1:nvars]   # interval length
    c = pf.scale_center                                         # center

    # Expand T_n((2*(x-c)/s)) to monomials in x
    function cheb_scaled_to_monomial(n, scale, center)
        tcoeffs = cheb_to_monomial_coeffs(n)  # T_n(X) -> X^k
        expanded = Dict{Int, Float64}()

        for (k, t) in enumerate(tcoeffs)
            k0 = k - 1
            factor = (2/scale)^k0  # scale inside expansion
            # expand (x-center)^k0 using binomial theorem
            for i in 0:k0
                expanded[i] = get(expanded, i, 0.0) + t * factor * binomial(k0,i) * (-center)^(k0 - i)
            end
        end
        return expanded
    end

    # Loop over all tensor-product Chebyshev terms
    for (coef, exps) in zip(pf.coeffs, pf.exponents)
        uni_expansions = [cheb_scaled_to_monomial(e, s[i], c[i]) for (i,e) in enumerate(exps)]
        iterators = [collect(pairs(v)) for v in uni_expansions]

        # Tensor-product over variables
        for prod_tuple in Iterators.product(iterators...)
            mono_exp = Tuple(k for (k, _) in prod_tuple)
            mono_coef = coef * prod(v for (_, v) in prod_tuple)
            monomial_coeffs[mono_exp] = get(monomial_coeffs, mono_exp, 0.0) + mono_coef
        end
    end

    return monomial_coeffs
end

pf = Polyfit([:x, :y], 6)

# generate sample data
data = []

for x in range(5, 8, length=12) .+ 0.2*im
    for y in range(3, 5, length=12) .+ 0.3*im
        fxy = x^6 + y^6 + 3x*y + 1e-5rand()
        push!(data, ([x,y], fxy))
    end
end

# fit the polynomial
fit!(pf, data)

monomials = polyfit_to_monomial(pf)

# show the fit in LaTeX (Chebyshev basis)
show(stdout, "text/latex", pf)