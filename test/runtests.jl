using JuliVirBootstrap
using Test

@testset "CFTData.jl" begin

    #ensure the relation between b and β does not change
    c1 = CentralCharge(:c, -1.1+.2im)
    b = c1.b
    c2 = CentralCharge(:b, b)
    @test c1.c == c2.c
    @test c1.β == c2.β

    #ensure the relation between p and P does not change
    left = 1
    right = 2
    V1 = Field(c1, :P, 0.5, diagonal=true)
    p = V1.P[left]
    V2 = Field(c1, :P, p, diagonal=true)
    @test V1.P == V2.P

    #ensure the keyword diagonal also works for fields given from Kac indices
    V1 = Field(c1, Kac=true, r=3, s=4, diagonal=true)
    @test V1.δ[left] == V1.δ[right]


    #ensure degenerate and diagonal work well together
    V1 = Field(c1, Kac=true, degenerate=true, r=2, s=5, diagonal=true)
    @test V1.δ[left] == V1.δ[right]

end

@testset "FourPointCorrelationFunctions" begin

    left=1
    right=2

    c = CentralCharge(:β, 1.2+.1*1im)
    V1 = Field(c, :Δ, big"0.23"+big".11"*im, diagonal=true)
    V2 = Field(c, :Δ, 3.43, diagonal=true)
    V3 = Field(c, :Δ, 0.13, diagonal=true)
    V4 = Field(c, :Δ, 1.3, diagonal=true)
    corr = FourPointCorrelation(c, V1, V2, V3, V4)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.Rmn(2, 1, corr, left),
                   0.31097697185245077-0.70523695127635733im, # value taken from Sylvain's code
                   atol=1e-8)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.computeCNmn(7, 2, 3, corr, left),
                   0.0019498393368877166+0.0026353877950837049im, # value taken from Sylvain's code
                   atol=1e-8)

end

@testset "FourPointBlocks" begin

    left=1;
    right=2;

    import JuliVirBootstrap.FourPointBlocksSphere.qfromx

c = CentralCharge(:b, (1.2+.1*1im)/sqrt(2))

q = JuliVirBootstrap.FourPointBlocksSphere.qfromx(0.05)

P = 0.23+.11im
P1 = 0.41+1.03im

V_chan = Field(c, :P, sqrt(2)*P, diagonal=true)
V_ext = Field(c, :P, P1/sqrt(2), diagonal=true)
VKac = Field(c, Kac=true, r=0, s=1//2, diagonal=true)

corr = FourPointCorrelation(c, [VKac, V_ext, VKac,VKac])
block = FourPointBlockSphere(corr, :s, V_chan)

h = evalpoly(16*q, JuliVirBootstrap.FourPointBlocksSphere.H_series(block, left))

@test isapprox(h, 0.9999955375834808 - 2.735498726466085e-6im, atol=1e-8) # value from Sylvain's code

setprecision(BigFloat, 128)

c = CentralCharge(:c, big"0.1")
V1 = Field(c, :Δ, 1, diagonal=true)
V2 = Field(c, :Δ, 2, diagonal=true)
V3 = Field(c, :Δ, 3, diagonal=true)
V4 = Field(c, :Δ, 4, diagonal=true)
corr = FourPointCorrelation(c, V1, V2, V3, V4)
V = Field(c, :Δ, big"0.5", diagonal=true)
block_s = FourPointBlockSphere(corr, :s, V, Nmax=50)
block_t = FourPointBlockSphere(corr, :t, V, Nmax=50)
block_u = FourPointBlockSphere(corr, :u, V, Nmax=50)
x=big"0.05"

# comparing to values from Sylvain's code
@test isapprox(block_chiral(x, block_s, left), big"1679.9121886897846270816517306779666391454311387606437056866150367", rtol = 1e-20)
@test isapprox(block_chiral(x, block_t, left), big"10841.2576587560092582414458316202779244541207",rtol = 1e-20)
@test isapprox(block_chiral(x, block_u, right), big"299.1846813850886027170806472436222922268361198327" -big"2026.4731585077561510727121083232012071890514123"*im, rtol = 1e-20)

setprecision(BigFloat, 64)
left = 1
right = 2

c = CentralCharge(:β, big".912" + .1im)
V1 = Field(c, Kac=true, r=1//2, s=0)
V2 = Field(c, Kac=true, r=3//2, s=2//3)

corr = FourPointCorrelation(c, [V1, V1, V2, V1])
block_s = FourPointBlockSphere(corr, :s, V1, Nmax=15)
block_t = FourPointBlockSphere(corr, :t, V1, Nmax=15)

z = 1e-8 + 1e-10im
Δ1 = V1.Δ[left]

@test abs(1-block_non_chiral(z, block_s)*z^Δ1*conj(z)^Δ1) < 1e-5
@test abs(1-block_non_chiral(1-z, block_t)*z^Δ1*conj(z)^Δ1) < 1e-5 # both blocks are close to one

setprecision(BigFloat, 256)

c = CentralCharge(:β, big(1.2 + .1im))
V1 = Field(c, Kac=true, r=1//2, s=0)
V2 = Field(c, Kac=true, r=3//2, s=2//3)

ϵ = big(1e-25)
V = Field(c, :Δ, 0.5, diagonal=true)
Vshiftedp = Field(c, :Δ, 0.5+ϵ, diagonal=true)
Vshiftedm = Field(c, :Δ, 0.5-ϵ, diagonal=true)

corr = FourPointCorrelation(c, [V1, V1, V2, V1])
block = FourPointBlockSphere(corr, :s, V, Nmax=50)
block_shiftedp = FourPointBlockSphere(corr, :s, Vshiftedp, Nmax=50)
block_shiftedm = FourPointBlockSphere(corr, :s, Vshiftedm, Nmax=50)

block_der = block_chiral(z, block, left, der=true)
block_der_manual = (block_chiral(z, block_shiftedp, left) - block_chiral(z, block_shiftedm, left))/(2*ϵ)

@test abs(block_der - block_der_manual) < 1e-6

c = CentralCharge(:β, big(.8 + .1im))
V1 = Field(c, Kac=true, r=1, s=1)
V2 = Field(c, Kac=true, r=1, s=1)
V3 = Field(c, Kac=true, r=0, s=1//2)
V4 = Field(c, Kac=true, r=0, s=3//2)
VΔ = Field(c, :Δ, 0.5, diagonal=true)

corr = FourPointCorrelation(c, [V1, V2, V3, V4])
corrΔ = FourPointCorrelation(c, [V1, V2, V3, VΔ])

ell = JuliVirBootstrap.FourPointBlocksSphere.ell(corr, 2, 1)
ellΔ = JuliVirBootstrap.FourPointBlocksSphere.ell(corrΔ, 2, 1)

# When all fields are degenerate
@test  isapprox(ell, 8.2808044631395529307 - 9.7096599503345083802im, rtol = 1e-8) # comparing with Sylvain's code
# When not all fields are degenerate
@test isapprox(ellΔ, 11.392850199938978801 - 7.6477614372039684265im, rtol = 1e-8) # comparing with Sylvain's code

c = CentralCharge(:β, big(1.2 + .1im))
V1 = Field(c, Kac=true, r=0, s=1)
V2 = Field(c, Kac=true, r=0, s=1//2)
V3 = Field(c, Kac=true, r=0, s=1)
V4 = Field(c, Kac=true, r=0, s=1//2)

V = Field(c, Kac=true, r=2, s=3)

x = 0.3 + 0.1im
Nmax = 26
corr = FourPointCorrelation(c, [V1, V2, V3, V4])
b(channel) = FourPointBlockSphere(corr, channel, V)
block_value(b) = block_non_chiral(x, b)

@test isapprox(block_value(b(:s)), -0.0062116451268237 + 0.0009314731786393im, rtol = 1e-5) # comparing with Sylvain's code
@test isapprox(block_value(b(:t)), -0.15830875034149818 - 0.130335270628475im, rtol = 1e-5) # comparing with Sylvain's code
@test isapprox(block_value(b(:u)), 296.0639291056886 - 16.68222738906im, rtol = 1e-3) # comparing with Sylvain's code:TODO precision problem

end

@testset "OnePointBlocks" begin
    left=1;
    right=2;

    import JuliVirBootstrap.FourPointBlocksSphere.qfromx
    c_torus = CentralCharge(:b, 1.2+.1*1im);
    c_sphere = CentralCharge(:b, (1.2+.1*1im)/sqrt(2))

    q = JuliVirBootstrap.FourPointBlocksSphere.qfromx(0.05)

    P = 0.23+.11im
    P1 = 0.41+1.03im
    V_torus_chan = Field(c_torus, :P, P, diagonal=true)
    δ_torus = V_torus_chan.δ[left]
    δ11_torus = Field(c_torus, Kac=true, r=1, s=1, diagonal=true).δ[left]
    V_torus_ext = Field(c_torus, :P, P1, diagonal=true)

    V_sphere_chan = Field(c_sphere, :P, sqrt(2)*P, diagonal=true)
    δ_sphere = V_sphere_chan.δ[left]
    δ21_sphere = Field(c_sphere, Kac=true, r=2, s=1, diagonal=true).δ[left]
    δ12_sphere = Field(c_sphere, Kac=true, r=1, s=2, diagonal=true).δ[left]
    V_sphere_ext = Field(c_sphere, :P, P1/sqrt(2), diagonal=true)
    VKac_sphere = Field(c_sphere, Kac=true, r=0, s=1//2, diagonal=true)

    corr_torus = OnePointCorrelation(c_torus, V_torus_ext)
    block_torus = OnePointBlockTorus(V_torus_chan)

    corr_sphere = FourPointCorrelation(c_sphere, [VKac_sphere, V_sphere_ext, VKac_sphere,VKac_sphere])
    block_sphere = FourPointBlockSphere(:s, V_sphere_chan)

    h1 = JuliVirBootstrap.OnePointBlocksTorus.H(q^2, 5, block_torus, corr_torus, left)
    h2 = JuliVirBootstrap.FourPointBlocksSphere.H(q, 5, block_sphere, corr_sphere, left)

    @test isapprox(h1, h2, atol=1e-12)
end
