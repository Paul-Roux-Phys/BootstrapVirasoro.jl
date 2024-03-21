#==================

SpecialFunctions.jl computes the special functions relevant for our applications in 2D CFT.

==================#

module SpecialFunctions

import SpecialFunctions: digamma as digamma_w_poles

export digamma

"""Regularised digamma function"""
function digamma(z)
    if real(z) > 0
        return digamma_w_poles(z)
    elseif imag(z) == 0 and real(z)%1 == 0
        return digamma_w_poles(1-z)
    else
        return digamma_w_poles(1-z) - π/tan(π*z)
    end
end

function log_double_gamma(beta, w)
end

function double_gamma(beta, w)
    return exp(log_double_gamma(beta, w))
end

end # end module
