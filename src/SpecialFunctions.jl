function jthetas(z, τ)
    t1 = zero(z)
    t2 = zero(z)
    t3 = zero(z)
    t4 = zero(z)

    Arblib.modular_theta!(t1, t2, t3, t4, z, τ)

    return (t1, t2, t3, t4)
end

jtheta1(z::Acb, τ::Acb) = jthetas(z, τ)[1]
jtheta2(z::Acb, τ::Acb) = jthetas(z, τ)[2]
jtheta3(z::Acb, τ::Acb) = jthetas(z, τ)[3]
jtheta4(z::Acb, τ::Acb) = jthetas(z, τ)[4]

jtheta1(τ::Acb) = jthetas(zero(τ), τ)[1]
jtheta2(τ::Acb) = jthetas(zero(τ), τ)[2]
jtheta3(τ::Acb) = jthetas(zero(τ), τ)[3]
jtheta4(τ::Acb) = jthetas(zero(τ), τ)[4]

function ellipticK(m::Acb)
    Arblib.elliptic_k!(zero(m), m)
end

function etaDedekind(τ::Acb)
    Arblib.modular_eta!(zero(τ), τ)
end

gamma(z::Acb) = Arblib.gamma!(zero(z), z)
gamma(z) = gamma(Acb(z))
digamma(z::Acb) = Arblib.digamma!(zero(z), z)
digamma(z) = digamma(Acb(z))

function digamma_reg(z::Acb)
    if real(z) > 0
        return digamma(z)
    elseif Arblib.contains_int(z)
        return digamma(1-z)
    else
        return digamma(1-z) - π*cot(π * z)
    end
end
