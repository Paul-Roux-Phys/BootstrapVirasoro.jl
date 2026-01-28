@testset "Normal residues" begin
    c = CC(c = 0.5)
    V1 = Field(c, r = 1, s = 0)
    V2 = Field(c, r = 2, s = 0)
    corr = Corr(V1, V2, V2, V1, 6)

    # litteral values are taken from Sylvain's code
    @test isapprox(corr.Rmn.left.s[3, 2], -3.11111e-7, atol = 1e-8)

    c = CC(c = big"0.1")
    V1 = Field(c, r = 1 // 2, s = 0)
    V2 = Field(c, r = 3 // 2, s = 0)
    corr = Corr(V1, V2, V2, V1, 4)

    @test isapprox(
        corr.CNmn.left.t[4, 1, 4],
        big"-0.00004010759186977254114462018639739857292228437487",
        atol = 1e-20,
    )
end

@testset "Regularised residues as limits of residues" begin
    import BootstrapVirasoro: computeRmns, Rmn_term, Rmn_zero_order, Rmn_term_vanishes
    setprecision(BigFloat, 256)
    c = CentralCharge(β = big"0.8" + big"0.1" * im)
    ϵ = 1 // big"10"^20

    d1 = CD(c, r = 0, s = 1)
    d2 = CD(c, r = 0, s = 1 // 2)
    d3 = CD(c, r = 2, s = 1 // 2)
    d3bar = CD(c, r = 2, s = -1 // 2)
    d4(ϵ) = CD(c, r = 2, s = 3 // 2 + ϵ)
    d4bar(ϵ) = CD(c, r = 2, s = -3 // 2 - ϵ)
    dls(ϵ) = (d1, d2, d3, d4(ϵ))
    drs(ϵ) = (d1, d2, d3bar, d4bar(ϵ))

    P3, P4 = d3.P, d4(ϵ).P
    P(r, s) = (r * c.β - s / c.β) / 2

    Δmax      = 12
    T         = Complex{BigFloat}
    DRs       = Matrix{T}(undef, (Δmax, Δmax))
    Pns       = Matrix{T}(undef, (Δmax, Δmax))
    factors   = Matrix{T}(undef, (Δmax, 2Δmax))

    Rl_reg    = computeRmns(Δmax, dls(0))[2]
    Rlϵ       = computeRmns(Δmax, dls(ϵ))[1]
    Rr_reg    = computeRmns(Δmax, drs(0))[2]
    Rrϵ       = computeRmns(Δmax, drs(ϵ))[1]

    Rlreg_1_2 = Rl_reg[1, 2]
    Rlϵ_1_2   = Rlϵ[1, 2] / (P3 - P4 + P(0, 1))

    @test isapprox(Rlreg_1_2, Rlϵ_1_2, rtol = 1e-18)

    Rreg_12 = Rl_reg[1, 12]
    barRreg_12 = Rr_reg[1, 12]
    Rϵ_12 = Rlϵ[1, 12]
    barRϵ_12 = Rrϵ[1, 12] # / (P3 - P4 + P(0, 1))

    Rratio_reg = (-1)^Rmn_zero_order(1, 12, dls(0)) * Rreg_12 / barRreg_12
    Rratio_ϵ = Rϵ_12 / barRϵ_12

    @test isapprox(Rratio_reg, Rratio_ϵ, rtol = 1e-18)
end

@testset "Zamolodchikov series" begin
    import BootstrapVirasoro: qfromx, eval_series

    c = CC(b = (1.2 + 0.1 * 1im) / sqrt(2))
    x = 0.05
    q = qfromx(x)
    P = 0.23 + 0.11im
    P1 = 0.41 + 1.03im

    V_chan = Field(c, P = sqrt(2) * P)
    V_ext = Field(c, P = P1 / sqrt(2))
    VKac = Field(c, r = 0, s = 1 // 2)

    corr = Correlation(VKac, V_ext, VKac, VKac, 12)
    b = CBlock(corr[:left], :s, V_chan[:left], 12)

    h = eval_series(b, complex(x))

    @test h ≈ Base.evalpoly(16q, b.coeffs)
    @test isapprox(h, 0.9999955375834808 - 2.735498726466085e-6im, atol = 1e-8) # value from Sylvain's code
end

c = CentralCharge(c = big"0.1")
Δmax = 40
V1 = CD(c, Δ = 1)
V2 = CD(c, Δ = 2)
V3 = CD(c, Δ = 3)
V4 = CD(c, Δ = 4)
co = Correlation(V1, V2, V3, V4, Δmax)
V = CD(c, Δ = big"0.5")
b_s = CBlock(co, :s, V, Δmax)
b_t = CBlock(co, :t, V, Δmax)
b_u = CBlock(co, :u, V, Δmax)
x = big"0.05" + big"0" * im

@testset "prefactors" begin
    cache = BootstrapVirasoro.PosCache(x, b_s)
    @test isapprox(
        cache.prefactor * (16cache.q)^b_s.chan_dim.δ,
        big"1813.32806084410414587456604",
        rtol = 1e-20,
    )
    cache = BootstrapVirasoro.PosCache(1 - x, b_t)
    @test isapprox(
        cache.prefactor * (16cache.q)^b_t.chan_dim.δ,
        big"0.07933043122650460823164",
        rtol = 1e-20,
    )
    cache = BootstrapVirasoro.PosCache(1 / big"1.5" + 0im, b_u)
    @test isapprox(
        cache.prefactor * (16cache.q)^b_u.chan_dim.δ,
        big"3.425385476422140172584631280130419",
        rtol = 1e-20,
    )
end

# comparing to values from Sylvain's code
@testset "block values" begin
    @test isapprox(b_s(complex(x)), big"1679.912188689784627081651", rtol = 1e-20)

    @test isapprox(
        b_t(complex(1 - x)),
        big"10841.257658755518924543654",
        rtol = 1e-20,
    )

    @test isapprox(
        b_u(1 / (big"1.5" + 0im)),
        big"67.6043205801146820843104",
        rtol = 1e-20,
    )
end

c = CC(β = big"0.8" + big"0.1" * im)
V = Field(c, r = 2, s = 3)
Δmax = 40

V1 = Field(c, r = 0, s = 1)
V2 = Field(c, r = 0, s = 1 // 2)
V3 = Field(c, r = 0, s = 1)
V4 = Field(c, r = 0, s = 1 // 2)

co = Correlation(V1, V2, V3, V4, Δmax)

x = big"0.3" + big"0.1" * im

@testset "Limit as z->0, z->1" begin
    cor = Correlation(V1, V1, V2, V1, 12)
    block_s = NCBlock(cor, :s, V1, 12)
    block_t = NCBlock(cor, :t, V1, 12)

    z = 1e-8 + 1e-10im
    Δ = V1[:left].Δ

    # both blocks should be close to one
    @test abs(1 - block_s(z) * z^Δ * conj(z)^Δ) < 1e-5
    @test abs(1 - block_t(z) * z^Δ * conj(z)^Δ) < 1e-5
end

@testset "Block derivatives" begin
    setprecision(BigFloat, 256)
    ϵ = 1e-25

    V0 = Field(c, P = big"0.5")
    Vp = Field(c, P = big"0.5" + ϵ)
    Vm = Field(c, P = big"0.5" - ϵ)

    b = CBlock(co[:left], :s, V0[:left], 40, true)
    b_p = CBlock(co[:left], :s, Vp[:left], 40)
    b_m = CBlock(co[:left], :s, Vm[:left], 40)

    function eval_series_der(x)
        series_der = BootstrapVirasoro.eval_series_der(b, x)
        byhand =
            (
                BootstrapVirasoro.eval_series(b_p, x) -
                BootstrapVirasoro.eval_series(b_m, x)
            ) / (2 * ϵ)
        series_der - byhand
    end

    @test abs(eval_series_der(big"0.3" + big"0.4" * im)) < big"1e-45"
    @test abs(eval_series_der(big"0.01")) < big"1e-45"
    @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-45"

    function eval_block_der(x)
        block_der = b(x, true)
        byhand = (b_p(x) - b_m(x)) / (2 * ϵ)
        block_der - byhand
    end

    @test abs(eval_block_der(big"0.3" + big"0.4" * im)) < big"1e-25"
    @test abs(eval_block_der(big"0.01" + 0im)) < big"1e-25"
    @test abs(eval_block_der(big"10" + big"0.01" * im)) < big"1e-25"
end

@testset "Logarithmic prefactor ell" begin
    import BootstrapVirasoro: ell

    l = BootstrapVirasoro.ell(co.fields, 2, 1)
    # comparing with Sylvain's code
    # When all fields are degenerate
    @test c.β * l ≈ 14.20389003630952076540 - 5.0517664348341790287im
end

@testset "Regularised blocks" begin
    b = CBlock(co[:left], :s, V[:left])
    ϵ = 1e-40
    dϵ = CD(c, δ = V[:left].δ + ϵ)
    dPϵ = CD(c, P = V[:left].P + ϵ)
    dminus = CD(c, r = V.r, s = -V[:left].s)
    bϵ = CBlock(co[:left], :s, dϵ)
    bPϵ = CBlock(co[:left], :s, dPϵ)
    bminus = CBlock(co[:left], :s, dminus)

    @test isapprox(
        bϵ(x),
        co.Rmn.left.s[V.r, V.s] / ϵ * bminus(x) + b(x),
        rtol = 1e-32,
    )

    @test isapprox(
        bPϵ(x),
        co.Rmn.left.s[V.r, V.s] / 2 / V[:left].P / ϵ * bminus(x) + b(x),
        rtol = 1e-33,
    )

    @test isapprox(
        b(x),
        big"0.51970140827959736684758007395822214" +
        big"0.5951179392484063703449815783272925" * im,
        rtol = 1e-25,
    )
end

@testset "Full logarithmic blocks" begin
    setprecision(BigFloat, 167)
    bl(channel) = NCBlock(co, channel, V, Δmax)

    # comparing with Sylvain's code
    @test isapprox(
        bl(:s)(x),
        big"-0.0004874448542139286383748521" -
        big"0.001382427546460296313978939" * im,
        rtol = 1e-20,
    )

    @test isapprox(
        bl(:t)(1 - x),
        big"-0.4855554411145733280520066" + big"0.1128101630322690069857833" * im,
        rtol = 1e-20,
    )

    @test isapprox(
        bl(:u)(1 / (5 - x)), # eval near 5 because near zero the numerical error
        # can be large
        big"-6.036027998137231362922e-6" + big"2.335826931375437289964e-5" * im,
        rtol = 1e-20,
    )
end

@testset "Accident. non-log from generic log" begin
    V1 = Field(c, r = 0, s = 1)
    V2 = Field(c, r = 0, s = 1 // 2)
    V3 = Field(c, r = 2, s = 1 // 2)
    V4 = Field(c, r = 2, s = 3 // 2)
    ϵ = big"1" // big"10"^20
    V4ϵ = Field(c, r = 2, s = 3 // 2 + ϵ)
    V = Field(c, r = 1, s = 4)
    Vop = Field(c, r = 1, s = -V.s)

    Δmax = 25

    co = Co(V1, V2, V3, V4, Δmax)
    coϵ = Co(V1, V2, V3, V4ϵ, Δmax)

    block(chan) = NCBlock(co, chan, V, Δmax)
    blockϵ(chan) = NCBlock(coϵ, chan, V, Δmax)

    bs, bsϵ, bt, btϵ = block(:s), blockϵ(:s), block(:t), blockϵ(:t)

    # @test isapprox(bs(x), bsϵ(x), rtol = 1e-18)
    @test isapprox(bt(x), btϵ(x), rtol = 1e-18)
end
