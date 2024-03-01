#==================

SpecialFunctions.jl computes the special functions relevant for our applications in 2D CFT.

==================#

using EllipticFunctions

function log_double_gamma(beta, w)
end

function double_gamma(beta, w)
    return exp(log_double_gamma(beta, w))
end
