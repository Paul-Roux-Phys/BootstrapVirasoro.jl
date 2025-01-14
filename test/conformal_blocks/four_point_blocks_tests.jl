@testset "Zamolodchikov series" begin
    import BootstrapVirasoro: qfromx, evaluate_series

    c = CentralCharge(:b, (1.2 + 0.1 * 1im) / sqrt(2))
    x = 0.05
    q = qfromx(x)
    P = 0.23 + 0.11im
    P1 = 0.41 + 1.03im

    V_chan = Field(c, :P, sqrt(2) * P, diagonal=true)
    V_ext = Field(c, :P, P1 / sqrt(2), diagonal=true)
    VKac = Field(c, r=0, s=1//2, diagonal=true)

    corr = Correlation(VKac, V_ext, VKac, VKac, 12)
    b = Block(corr, :s, V_chan, :left, 12)

    h = evaluate_series(b, 16q)

    @test isapprox(h, 0.9999955375834808 - 2.735498726466085e-6im, atol=1e-8) # value from Sylvain's code

end

@testset "Chiral blocks" begin
    c = CentralCharge(:c, big"0.1")
    V1 = Field(c, :Δ, 1, diagonal=true)
    V2 = Field(c, :Δ, 2, diagonal=true)
    V3 = Field(c, :Δ, 3, diagonal=true)
    V4 = Field(c, :Δ, 4, diagonal=true)
    co = Correlation(V1, V2, V3, V4, 40)
    V = Field(c, :Δ, big"0.5", diagonal=true)
    b_s = Block(co, :s, V, :left, 40)
    b_t = Block(co, :t, V, :left, 40)
    b_u = Block(co, :u, V, :left, 40)
    x=big"0.05"

    # comparing to values from Sylvain's code
    @test isapprox(
        evaluate(b_s, x),
        big"1679.9121886897846270816517306779666391454311387606437056866150367",
        rtol = 1e-20
    )

    @test isapprox(
        evaluate(b_t, x),
        big"10841.257658755518924543654668282368582",
        rtol = 1e-20
    )

    @test isapprox(
        evaluate(b_u, x),
        big"299.1843849134281379909218349390129720213727" -
            big"2026.47316197194594006925827438887723275377012"*im,
        rtol = 1e-20
    )
end


@testset "Non Chiral Blocks" begin

    c = CentralCharge(:β, big"0.8"+big"0.1"*im)
    V = Field(c, r=2, s=3)
    Nmax = 26

    V1 = Field(c, r=0, s=1)
    V2 = Field(c, r=0, s=1//2)
    V3 = Field(c, r=0, s=1)
    V4 = Field(c, r=0, s=1//2)
    VΔ = Field(c, :Δ, big"0.5", diagonal=true)

    co = Correlation(V1, V2, V3, V4, Nmax)
    coΔ = Correlation(V1, V2, V3, VΔ, Nmax)

    x = big"0.3"+big"0.1"*im
    
    @testset "Limit as z->0, z->1" begin

        cor = Correlation(V1, V1, V2, V1, 12)
        block_s = Block(cor, :s, V1, 12)
        block_t = Block(cor, :t, V1, 12)

        z = 1e-8 + 1e-10im
        Δ = V1.Δ[:left]

        # both blocks should be close to one
        @test abs(1 - evaluate(block_s, z) * z^Δ * conj(z)^Δ) < 1e-5
        @test abs(1 - evaluate(block_t, 1 - z) * z^Δ * conj(z)^Δ) < 1e-5
    end

    @testset "Block derivatives" begin
        setprecision(BigFloat, 256)
        ϵ = 1e-25

        V0 = Field(c, :P, big"0.5", diagonal=true)
        Vp = Field(c, :P, big"0.5" + ϵ, diagonal=true)
        Vm = Field(c, :P, big"0.5" - ϵ, diagonal=true)

        b = Block(co, :s, V0, :left, 40, der=true)
        b_p = Block(co, :s, Vp, :left, 40)
        b_m = Block(co, :s, Vm, :left, 40)

        function eval_series_der(x)
            q = 16qfromx(x)
            series_der = evaluate_series(b, q, der=true)
            byhand = (evaluate_series(b_p, q) - evaluate_series(b_m, q)) / (2 * ϵ)
            series_der - byhand
        end

        @test abs(eval_series_der(big"0.3" + big"0.4" * im)) < big"1e-45"
        @test abs(eval_series_der(big"0.01")) < big"1e-45"
        @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-45"

        function eval_block_der(x)
            block_der = evaluate(b, x, der=true)
            byhand = (evaluate(b_p, x) - evaluate(b_m, x)) / (2 * ϵ)
            block_der - byhand
        end

        @test abs(eval_block_der(big"0.3" + big"0.4" * im)) < big"1e-25"
        @test abs(eval_block_der(big"0.01")) < big"1e-25"
        @test abs(eval_block_der(big"10" + big"0.01" * im)) < big"1e-25"
    end


    @testset "Logarithmic prefactor ell" begin
        import BootstrapVirasoro: ell

        l = ell(co.fields, :s, 2, 1)
        lΔ = ell(coΔ.fields, :s, 2, 1)

        # comparing with Sylvain's code
        # When all fields are degenerate
        @test isapprox(l, 14.20389003630952076540 - 5.0517664348341790287im, rtol=1e-15)
        # When not all fields are degenerate
        @test isapprox(lΔ, 9.442859099125287026870 - 8.0848336893143160748im, rtol=1e-15)
    end

    @testset "Regularised blocks" begin
        b = Block(co, :s, V, :left)
        V12 = Field(c, r=1, s=2)
        b2 = Block(coΔ, :t, V12, :left, Nmax)

        @test isapprox(
            evaluate(b, x),
            big"0.51970140827959736684758007395822214" + 
                big"0.5951179392484063703449815783272925"*im,
            rtol=1e-25
        )

        @test isapprox(
            evaluate(b2, x),
            big"1.164057039115389253869277410981757" + 
                big"-0.07162461483111554418789576414972262"*im,
            rtol=1e-25
        )
    end

    @testset "Full logarithmic blocks" begin
        bl(channel) = Block(co, channel, V, Nmax)

        # comparing with Sylvain's code
        @test isapprox(
            evaluate(bl(:s), x),
            big"-0.0004874448542139286383748521" - big"0.001382427546460296313978939"*im,
            rtol = 1e-20
        )
        @test isapprox(
            evaluate(bl(:t), x),
            big"-0.4855554411145733280520066" + big"0.1128101630322690069857833"*im,
            rtol = 1e-20
        )
        @test isapprox(
            evaluate(bl(:u), 5-x), # evaluate near 5 because near zero the numerical error
                                   # can be large
            big"-6.036027998137231362922e-6" + big"2.335826931375437289964e-5"*im,
            rtol = 1e-20
        )
    end

    @testset "Accident. non-log from generic log" begin
        V1 = Field(c, r=0, s=1)
        V2 = Field(c, r=0, s=1//2)
        V3 = Field(c, r=2, s=1//2)
        V_4(ϵ) = Field(c, r=2, s=3//2 + ϵ)

        V = Field(c, r=1, s=12)

        corr(ϵ) = Correlation(V1, V2, V3, V_4(ϵ), Nmax)
        ϵ = big"1" // big"10"^20

        block(chan, ϵ) = Block(corr(ϵ), chan, V, Nmax)

        redirect_stderr(devnull) do
            @test isapprox(
                evaluate(block(:s, 0), x),
                evaluate(block(:s, ϵ), x),
                rtol = 1e-18
            )

            @test isapprox(
                evaluate(block(:t, 0), x),
                evaluate(block(:t, ϵ), x),
                rtol = 1e-18
            )
        end
    end
end
