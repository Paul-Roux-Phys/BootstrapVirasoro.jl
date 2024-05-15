#==================

SpecialFunctions.jl computes the special functions relevant for our applications in 2D CFT.

==================#

module SpecialFunctions

using SpecialFunctions, ArbNumerics, Memoization
export Barnes_G, log_double_Gamma, double_Gamma

# the SpecialFunctions package has no arbitrary-precision complex-variable gamma function, however the ArbNumerics does. We use this, and convert to a Complex{BigFloat}

function log_Γ(z)
    return Complex{BigFloat}(lgamma(ArbComplex(z)))
end

function Γ(z)
    return Complex{BigFloat}(gamma(ArbComplex(z)))
end

function ψ(z)
    return Complex{BigFloat}(digamma(ArbComplex(z)))
end

function trigamma(z)
    return Complex{BigFloat}(polygamma(ArbComplex(1), ArbComplex(z)))
end

function polyΓ(n, z)
    return Complex{BigFloat}(polygamma(ArbComplex(n), ArbComplex(z)))
end
using QuadGK # numerical integration
using Symbolics, Memoization

export digamma_reg

"""Regularised digamma function"""
function digamma_reg(z)
    if real(z) > 0
        return ψ(z)
    elseif imag(z) == 0 && real(z)%1 == 0
        return ψ(1-z)
    else
        return ψ(1-z) - big(π)/tan(π*z)
    end
end

function integrand_C(x, τ)
    x = big(x)
    return exp((1-τ)*x)/(2*sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)*(exp(x)/(2*sinh(x))+1-τ/2)
end

function modular_C(τ)
    P = precision(BigFloat, base=10)
    #temporarily increase precision to avoid artificial divergence around zero
    setprecision(BigFloat, base=10, Int(floor(1.3*P)))
    cutoff = big(10^(-P/5)) # to prevent artificial divergence around zero
    tol = big(10)^P
    value, error = quadgk(x -> integrand_C(x, τ), cutoff, big(Inf), rtol = tol, order=21)
    C0 = (2/τ - 3//2 + τ/6)*cutoff + (5//12 - 1/τ + τ/12)*cutoff^2 + (4/(9*τ) - 2//9 + 1//54*τ - 1//270*τ^3)*cutoff^3
    setprecision(BigFloat, base=10, P)
    return 1/(2*τ)*log(2*big(π)) - value - C0
end

function integrand_D(x, τ)
    x = big(x)
    return x*exp((1-τ)*x)/(sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)
end

function modular_D(τ)
    P = precision(BigFloat, base=10)
    #temporarily increase precision to avoid artificial divergence around zero
    setprecision(BigFloat, base=10, Int(floor(1.3*P)))
    cutoff = big(10^(-P/5)) # to prevent artificial divergence around zero
    tol = big(10)^P
    value, error = quadgk( x -> integrand_D(x, τ), big(0), big(Inf), rtol = tol, order=21)
    setprecision(BigFloat, base=10, P)
    return value
end

@memoize function modular_coeff_a(τ)
    return 1/2*τ*log(big(2)*π*τ) + 1/2*log(τ) - τ*modular_C(τ)
end

@memoize function modular_coeff_b(τ)
    return -τ*log(τ) - τ^2*modular_D(τ)
end

function log_Barnes_GN(N, z, τ)
    term1 = - log(τ) - log_Γ(z)
    term2 = modular_coeff_a(τ)*z/τ + modular_coeff_b(τ)*z^2/(2*τ^2)
    term3 = sum(log_Γ(m*τ) - log_Γ(z+m*τ) + z*ψ(m*τ)+z^2/2*trigamma(m*τ) for m in 1:N)
    return term1 + term2 + term3
end

function polynomial_Pn(n, z, τ)
    if n == 1
        return 1//6
    else
        term1 = z^(n-1)/factorial(big(n+2))
        summand(k) = ((1+τ)^(k+2) - 1 - τ^(k+2))/(factorial(big(k+2))*τ) * polynomial_Pn(n-k, z, τ)
        return term1 - sum(summand(k) for k in 1:n-1)
    end
end

function rest_RMN(M, N, z, τ)
    return sum(factorial(big(k-1))*(-τ)^(-k-1)*polynomial_Pn(k, z, -τ)/N^k for k in 1:M)
end

"""Numerical approximation of the logarithm of Barne's G-function, up to a given tolerance"""
function log_Barnes_G(z, τ, tol)
    z = complex(z)
    d = -log(tol)/log(10)
    M = BigInt(floor(0.7*log(10)/log(20)*d))
    N = 20*M
    return log_Barnes_GN(N, z, τ) + z^3*rest_RMN(M, N, z, τ)
end

function Barnes_G(z, τ, tol)
    return exp(log_Barnes_G(z, τ, tol))
end

function log_Gamma_2(w, β, tol)
    β = real(β-1/β) < 0 ? 1/β : β # change β -> 1/β if needed
    return w/(2*β)*log(big(2)*π) + (w/2*(w-β-1/β)+1)*log(β) - log_Barnes_G(w/β, 1/β^2, tol)
end

"""
        log_double_Gamma(w, β, tol)

Compute the logarithm of the double gamma function Γ_β(w, β) with precision tol

"""
function log_double_Gamma(w, β, tol)
    return log_Gamma_2(w, β, tol) - log_Gamma_2((β+1/β)/2, β, tol)
end

"""
        double_Gamma(w, β, tol)

Compute the double gamma function Γ_β(w, β) with precision tol

"""
function double_Gamma(w, β, tol)
    exp(log_double_Gamma(w, β, tol))
end

end # end module
