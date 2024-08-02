@testset "Zamolodchikov series" begin
    import JuliVirBootstrap.ConformalBlocks: qfromx

    c = CentralCharge(:b, (1.2 + 0.1 * 1im) / sqrt(2))
    x = 0.05
    q = qfromx(x)
    P = 0.23 + 0.11im
    P1 = 0.41 + 1.03im

    V_chan = Field(c, :P, sqrt(2) * P, diagonal=true)
    V_ext = Field(c, :P, P1 / sqrt(2), diagonal=true)
    VKac = Field(c, Kac=true, r=0, s=1//2, diagonal=true)

    corr = FourPointCorrelation(c, VKac, V_ext, VKac, VKac, 20)
    b = FourPointBlock(20, c, corr, :s, V_chan)

    h = evaluate_series(b, 16q, left)

    @test isapprox(h, 0.9999955375834808 - 2.735498726466085e-6im, atol=1e-8) # value from Sylvain's code

end

@testset "Chiral blocks" begin
    c = CentralCharge(:c, big"0.1")
    V1 = Field(c, :Δ, 1, diagonal=true)
    V2 = Field(c, :Δ, 2, diagonal=true)
    V3 = Field(c, :Δ, 3, diagonal=true)
    V4 = Field(c, :Δ, 4, diagonal=true)
    corr = FourPointCorrelation(c, V1, V2, V3, V4, 40)
    V = Field(c, :Δ, big"0.5", diagonal=true)
    block_s = FourPointBlock(40, c, corr, :s, V)
    block_t = FourPointBlock(40, c, corr, :t, V)
    block_u = FourPointBlock(40, c, corr, :u, V)
    x=big"0.05"

    # comparing to values from Sylvain's code
    @test isapprox(
        evaluate_chiral(c, corr, block_s, x, left),
        big"1679.9121886897846270816517306779666391454311387606437056866150367",
        rtol = 1e-20
    )

    @test isapprox(
        evaluate_chiral(c, corr, block_t, x, left),
        big"10841.257658755518924543654668282368582",
        rtol = 1e-20
    )

    @test isapprox(
        evaluate_chiral(c, corr, block_u, x, right),
        big"299.1843849134281379909218349390129720213727" -
            big"2026.47316197194594006925827438887723275377012"*im,
        rtol = 1e-20
    )
end

@testset "Non Chiral Blocks" begin
    c = CentralCharge(:β, big".912" + .1im)
    V1 = Field(c, Kac=true, r=1//2, s=0)
    V2 = Field(c, Kac=true, r=3//2, s=2//3)

    corr = FourPointCorrelation(c, V1, V1, V2, V1, 20)
    block_s = FourPointBlock(15, c, corr, :s, V1)
    block_t = FourPointBlock(15, c, corr, :t, V1)

    z = 1e-8 + 1e-10im
    Δ1 = V1.Δ[left]

    # both blocks should be close to one
    # @test abs(1-evaluate_non_chiral(z, block_s)*z^Δ1*conj(z)^Δ1) < 1e-5
    # @test abs(1-evaluate_non_chiral(1-z, block_t)*z^Δ1*conj(z)^Δ1) < 1e-5 
end

@testset "Block derivatives" begin
    # setprecision(BigFloat, 256)

    # c = CentralCharge(:β, big(1.2 + .1im))
    # V1 = Field(c, Kac=true, r=1//2, s=0)
    # V2 = Field(c, Kac=true, r=3//2, s=2//3)

    # ϵ = big(1e-25)
    # V = Field(c, :Δ, 0.5, diagonal=true)
    # Vshiftedp = Field(c, :Δ, 0.5+ϵ, diagonal=true)
    # Vshiftedm = Field(c, :Δ, 0.5-ϵ, diagonal=true)

    # corr = FourPointCorrelation(c, [V1, V1, V2, V1])
    # block = FourPointBlockSphere(corr, :s, V, Nmax=50)
    # block_shiftedp = FourPointBlockSphere(corr, :s, Vshiftedp, Nmax=50)
    # block_shiftedm = FourPointBlockSphere(corr, :s, Vshiftedm, Nmax=50)

    # block_der = block_chiral(z, block, left, der=true)
    # block_der_manual = (block_chiral(z, block_shiftedp, left) - block_chiral(z, block_shiftedm, left))/(2*ϵ)

    # @test abs(block_der - block_der_manual) < 1e-6
end

@testset "Logarithmic blocks" begin
    # c = CentralCharge(:β, big(.8 + .1im))
    # V1 = Field(c, Kac=true, r=1, s=1)
    # V2 = Field(c, Kac=true, r=1, s=1)
    # V3 = Field(c, Kac=true, r=0, s=1//2)
    # V4 = Field(c, Kac=true, r=0, s=3//2)
    # VΔ = Field(c, :Δ, 0.5, diagonal=true)

    # corr = FourPointCorrelation(c, [V1, V2, V3, V4])
    # corrΔ = FourPointCorrelation(c, [V1, V2, V3, VΔ])

    # ell = JuliVirBootstrap.FourPointBlocksSphere.ell(corr, 2, 1)
    # ellΔ = JuliVirBootstrap.FourPointBlocksSphere.ell(corrΔ, 2, 1)

    # # When all fields are degenerate
    # @test  isapprox(ell, 8.2808044631395529307 - 9.7096599503345083802im, rtol = 1e-8) # comparing with Sylvain's code
    # # When not all fields are degenerate
    # @test isapprox(ellΔ, 11.392850199938978801 - 7.6477614372039684265im, rtol = 1e-8) # comparing with Sylvain's code

    # c = CentralCharge(:β, big(1.2 + .1im))
    # V1 = Field(c, Kac=true, r=0, s=1)
    # V2 = Field(c, Kac=true, r=0, s=1//2)
    # V3 = Field(c, Kac=true, r=0, s=1)
    # V4 = Field(c, Kac=true, r=0, s=1//2)

    # V = Field(c, Kac=true, r=2, s=3)

    # x = 0.3 + 0.1im
    # Nmax = 26
    # corr = FourPointCorrelation(c, [V1, V2, V3, V4])
    # b(channel) = FourPointBlockSphere(corr, channel, V)
    # block_value(b) = block_non_chiral(x, b)

    # @test isapprox(block_value(b(:s)), -0.0062116451268237 + 0.0009314731786393im, rtol = 1e-5) # comparing with Sylvain's code
    # @test isapprox(block_value(b(:t)), -0.15830875034149818 - 0.130335270628475im, rtol = 1e-5) # comparing with Sylvain's code
    # @test isapprox(block_value(b(:u)), 296.0639291056886 - 16.68222738906im, rtol = 1e-3) # comparing with Sylvain's code:TODO precision problem

end