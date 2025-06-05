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

    h = evaluate_series(b, complex(x))

    @test h ≈ Base.evalpoly(16q, b._coefficients)
    @test isapprox(h, 0.9999955375834808 - 2.735498726466085e-6im, atol=1e-8) # value from Sylvain's code
end

c = CentralCharge(:c, big"0.1")
Nmax = 40
V1 = ConformalDimension(c, :Δ, 1)
V2 = ConformalDimension(c, :Δ, 2)
V3 = ConformalDimension(c, :Δ, 3)
V4 = ConformalDimension(c, :Δ, 4)
co = Correlation(V1, V2, V3, V4, Nmax)
V = ConformalDimension(c, :Δ, big"0.5")
b_s = Block(co, :s, V, Nmax)
b_t = Block(co, :t, V, Nmax)
b_u = Block(co, :u, V, Nmax)
x = big"0.05" + big"0" * im

cache = BootstrapVirasoro.PositionCache(1 / big"1.5" + 0im, b_u)
cache.prefactor * (16cache.q)^b_u.channel_dimension.δ

@testset "prefactors" begin
    import BootstrapVirasoro: PositionCache
    cache = PositionCache(x, b_s)
    @test isapprox(
        cache.prefactor * (16cache.q)^b_s.channel_dimension.δ,
        big"1813.32806084410414587456604",
        rtol=1e-20
    )
    cache = PositionCache(1 - x, b_t)
    @test isapprox(
        cache.prefactor * (16cache.q)^b_t.channel_dimension.δ,
        big"0.07933043122650460823164",
        rtol=1e-20
    )
    cache = BootstrapVirasoro.PositionCache(1 / big"1.5" + 0im, b_u)
    @test isapprox(
        cache.prefactor * (16cache.q)^b_u.channel_dimension.δ,
        big"3.425385476422140172584631280130419",
        rtol=1e-20
    )
end

# comparing to values from Sylvain's code
@testset "block values" begin
    @test isapprox(
        evaluate(b_s, complex(x)),
        big"1679.912188689784627081651",
        rtol=1e-20
    )

    @test isapprox(
        evaluate(b_t, complex(1 - x)),
        big"10841.257658755518924543654",
        rtol=1e-20
    )

    @test isapprox(
        evaluate(b_u, 1 / (big"1.5" + 0im)),
        big"67.6043205801146820843104",
        rtol=1e-20
    )
end

c = CentralCharge(:β, big"0.8" + big"0.1" * im)
V = Field(c, r=2, s=3)
Nmax = 40

V1 = Field(c, r=0, s=1)
V2 = Field(c, r=0, s=1 // 2)
V3 = Field(c, r=0, s=1)
V4 = Field(c, r=0, s=1 // 2)

co = Correlation(V1, V2, V3, V4, Nmax)

x = big"0.3" + big"0.1" * im

@testset "Limit as z->0, z->1" begin
    cor = Correlation(V1, V1, V2, V1, 12)
    block_s = Block(cor, :s, V1, 12)
    block_t = Block(cor, :t, V1, 12)

    z = 1e-8 + 1e-10im
    Δ = V1.Δ[:left]

    # both blocks should be close to one
    @test abs(1 - evaluate(block_s, z) * z^Δ * conj(z)^Δ) < 1e-5
    @test abs(1 - evaluate(block_t, z) * z^Δ * conj(z)^Δ) < 1e-5
end

@testset "Block derivatives" begin
    import BootstrapVirasoro: evaluate_series, evaluate_series_der
    setprecision(BigFloat, 256)
    ϵ = 1e-25

    V0 = Field(c, :P, big"0.5", diagonal=true)
    Vp = Field(c, :P, big"0.5" + ϵ, diagonal=true)
    Vm = Field(c, :P, big"0.5" - ϵ, diagonal=true)

    b = Block(co, :s, V0, :left, 40, der=true)
    b_p = Block(co, :s, Vp, :left, 40)
    b_m = Block(co, :s, Vm, :left, 40)

    function eval_series_der(x)
        series_der = evaluate_series_der(b, x)
        byhand = (evaluate_series(b_p, x) - evaluate_series(b_m, x)) / (2 * ϵ)
        series_der - byhand
    end

    @test abs(eval_series_der(big"0.3" + big"0.4" * im)) < big"1e-45"
    @test abs(eval_series_der(big"0.01")) < big"1e-45"
    @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-45"

    function eval_block_der(x)
        block_der = evaluate_der(b, x)
        byhand = (evaluate(b_p, x) - evaluate(b_m, x)) / (2 * ϵ)
        block_der - byhand
    end

    @test abs(eval_block_der(big"0.3" + big"0.4" * im)) < big"1e-25"
    @test abs(eval_block_der(big"0.01" + 0im)) < big"1e-25"
    @test abs(eval_block_der(big"10" + big"0.01" * im)) < big"1e-25"
end

@testset "Logarithmic prefactor ell" begin
    import BootstrapVirasoro: ell

    l = ell(co.fields, :s, 2, 1)
    # comparing with Sylvain's code
    # When all fields are degenerate
    @test c.β * l ≈ 14.20389003630952076540 - 5.0517664348341790287im
end

@testset "Regularised blocks" begin
    b = Block(co, :s, V, :left)
    ϵ = 1e-40
    dϵ = ConformalDimension(c, δ=V.δ[:left] + ϵ)
    dPϵ = ConformalDimension(c, P=V.P[:left] + ϵ)
    dminus = ConformalDimension(c, r=V.r, s=-V.dims[:left].s)
    bϵ = Block(co[:left], :s, dϵ)
    bPϵ = Block(co[:left], :s, dPϵ)
    bminus = Block(co[:left], :s, dminus)

    @test isapprox(
        evaluate(bϵ, x),
        co._Rmn[:left][:s][V.r, V.s] / ϵ * evaluate(bminus, x) +
        evaluate(b, x),
        rtol=1e-32
    )

    @test isapprox(
        evaluate(bPϵ, x),
        co._Rmn[:left][:s][V.r, V.s] / 2 / V.P[:left] / ϵ * evaluate(bminus, x) +
        evaluate(b, x),
        rtol=1e-33
    )

    V12 = Field(c, r=1, s=2)
    @test isapprox(
        evaluate(b, x),
        big"0.51970140827959736684758007395822214" +
        big"0.5951179392484063703449815783272925" * im,
        rtol=1e-25
    )
end

@testset "Full logarithmic blocks" begin
    bl(channel) = Block(co, channel, V, Nmax)

    evaluate(bl(:s), x)

    # comparing with Sylvain's code
    @test isapprox(
        evaluate(bl(:s), x),
        big"-0.0004874448542139286383748521" - big"0.001382427546460296313978939" * im,
        rtol=1e-20
    )

    @test isapprox(
        evaluate(bl(:t), 1 - x),
        big"-0.4855554411145733280520066" + big"0.1128101630322690069857833" * im,
        rtol=1e-20
    )

    @test isapprox(
        evaluate(bl(:u), 1 / (5 - x)), # evaluate near 5 because near zero the numerical error
        # can be large
        big"-6.036027998137231362922e-6" + big"2.335826931375437289964e-5" * im,
        rtol=1e-20
    )
end

@testset "Accident. non-log from generic log" begin
    V1 = Field(c, r=0, s=1)
    V2 = Field(c, r=0, s=1 // 2)
    V3 = Field(c, r=2, s=1 // 2)
    V_4 = Field(c, r=2, s=3 // 2)
    ϵ = big"1" // big"10"^20
    V_4ϵ = Field(c, r=2, s=3 // 2 + ϵ)
    V = Field(c, r=1, s=12)
    Vop = Field(c, r=1, s=-V.s)

    correl = Correlation(V1, V2, V3, V_4, Nmax)
    correlϵ = Correlation(V1, V2, V3, V_4ϵ, Nmax)
    block(chan) = Block(correl, chan, V, Nmax)
    blockϵ(chan) = Block(correlϵ, chan, V, Nmax)

    bs, bsϵ, bt, btϵ = redirect_stderr(devnull) do
        block(:s), blockϵ(:s), block(:t), blockϵ(:t)
    end

    @test isapprox(bs(x), bsϵ(x), rtol=1e-18)
    @test isapprox(bt(x), btϵ(x), rtol=1e-18)
end

@testset "Interchiral" begin
    V = Field(c, r=2, s=1 // 2)
    b = Block(co, :s, V, interchiral=true, Δmax=30)
    @test b.shifts[1] == 1
    @test isapprox(
        b.shifts[3], big"5.464439142274867e-17" + big"1.7406391176219884e-17" * im,
        rtol=1e-15
    )
end
