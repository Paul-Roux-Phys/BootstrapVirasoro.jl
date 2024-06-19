#==================

SpecialFunctions.jl computes the special functions relevant for our applications in 2D CFT.

==================#

module SpecialFunctions

import SpecialFunctions as SF
using Memoization
using ArbNumerics # the SpecialFunctions package has no arbitrary-precision complex-variable gamma function, however the ArbNumerics does. We use this, and convert to a Complex{BigFloat}
using QuadGK # numerical integration

export loggamma, gamma
export digamma_reg, Barnes_G, logdoublegamma, doublegamma

"""
    cotpi(x) = cot(π * x)
"""
cotpi(x) = SF._cotpi(x)

ArbNumerics.loggamma = ArbNumerics.lgamma # rename ArbNumerics' loggamma function

for f in (:gamma, :digamma)
    @eval $f(z::Union{Real, Complex{Float64}}) = SF.$f(z)
    @eval $f(z::Complex{BigFloat}) = Complex{BigFloat}(ArbNumerics.$f(ArbComplex(z, bits=precision(BigFloat))))
end

trigamma(z::Union{Float64, Complex{Float64}}) = SF.trigamma(z)
trigamma(z::Union{BigFloat, Complex{BigFloat}}) = BigFloat(ArbNumerics.polygamma(ArbComplex(1, bits=precision(BigFloat)), ArbComplex(z, bits=precision(BigFloat))))

loggamma(z::Union{Real, Complex{Float64}}) = SF.loggamma(z)
loggamma(z::Complex{BigFloat}) = Complex{BigFloat}(ArbNumerics.lgamma(ArbComplex(z, bits=precision(BigFloat))))


polygamma(n, z::Union{Real, Complex{Float64}}) = SF.polygamma(n, z)
polygamma(n, z::Complex{BigFloat}) = Complex{BigFloat}(polygamma(ArbComplex(n), ArbComplex(z, bits=precision(BigFloat))))

"""Regularised digamma function"""
function digamma_reg(z)
    if real(z) > 0
        return digamma(z)
    elseif isreal(z) && real(z) < 0 && real(z)%1 == 0
        return digamma(1-z)
    else
        return digamma(1-z) - oftype(z, π)*cotpi(z)
    end
end

function integrandC(x, τ)
    x = big(x)
    return exp((1-τ)*x)/(2*sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)*(exp(x)/(2*sinh(x))+1-τ/2)
end

function modularC(τ)
    P = precision(BigFloat, base=10)
    #temporarily increase precision to avoid artificial divergence around zero
    setprecision(BigFloat, base=10, Int(floor(1.3*P)))
    cutoff = big(10^(-P/5)) # to prevent artificial divergence around zero
    tol = big(10)^P
    value, error = quadgk(x -> integrandC(x, τ), cutoff, big(Inf), rtol = tol, order=21)
    C0 = (2/τ - 3//2 + τ/6)*cutoff + (5//12 - 1/τ + τ/12)*cutoff^2 + (4/(9*τ) - 2//9 + 1//54*τ - 1//270*τ^3)*cutoff^3
    setprecision(BigFloat, base=10, P)
    return 1/(2*τ)*log(2*oftype(τ, π)) - value - C0
end

function integrandD(x, τ)
    x = big(x)
    return x*exp((1-τ)*x)/(sinh(x)*sinh(τ*x)) - exp(-2*x)/(τ*x)
end

function modularD(τ)
    P = precision(BigFloat, base=10)
    #temporarily increase precision to avoid artificial divergence around zero
    setprecision(BigFloat, base=10, Int(floor(1.3*P)))
    cutoff = big(10^(-P/5)) # to prevent artificial divergence around zero
    tol = big(10)^P
    value, error = quadgk( x -> integrandD(x, τ), big(0), big(Inf), rtol = tol, order=21)
    setprecision(BigFloat, base=10, P)
    return value
end

@memoize function modularcoeff_a(τ)
    return 1/2*τ*log(2*oftype(τ, π)*τ) + 1/2*log(τ) - τ*modularC(τ)
end

@memoize function modularcoeff_b(τ)
    return -τ*log(τ) - τ^2*modularD(τ)
end

function log_Barnes_GN(N, z, τ)
    res = 0
    res += - log(τ) - loggamma(z)
    res += modularcoeff_a(τ)*z/τ + modularcoeff_b(τ)*z^2/(2*τ^2)
    res += sum(loggamma(m*τ) - loggamma(z+m*τ) + z*digamma(m*τ)+z^2/2*trigamma(m*τ) for m in 1:N)
    return res
end

@memoize function factorial_big(n)::BigInt
    return factorial(big(n))
end

@memoize function polynomial_Pn(n, z, τ)
    if n == 1
        return 1//6
    else
        term1 = z^(n-1)/factorial_big(n+2)
        summand(k) = ((1+τ)^(k+2) - 1 - τ^(k+2))/(factorial_big(k+2)*τ) * polynomial_Pn(n-k, z, τ)
        return term1 - sum(summand(k) for k in 1:n-1)
    end
end

function rest_RMN(M, N, z, τ)
    return sum(factorial_big(k-1)*(-τ)^(-k-1)*polynomial_Pn(k, z, -τ)/N^k for k in 1:M)
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
    return w/(2*β)*log(2*oftype(β, π)) + (w/2*(w-β-1/β)+1)*log(β) - log_Barnes_G(w/β, 1/β^2, tol)
end

"""
        logdoublegamma(w, β, tol) = Γ_β(w)

Compute the logarithm of the double gamma function Γ_β(w, β) with precision tol
"""
function logdoublegamma(w, β, tol)
    return log_Gamma_2(w, β, tol) - log_Gamma_2((β+1/β)/2, β, tol)
end

"""
        doublegamma(w, β, tol)

Compute the double gamma function Γ_β(w, β) with precision tol

"""
function doublegamma(w, β, tol)
    exp(logdoublegamma(w, β, tol))
end

end # end module
