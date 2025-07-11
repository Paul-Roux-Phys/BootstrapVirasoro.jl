using BootstrapVirasoro: islogarithmic, PositionCache

setprecision(BigFloat, 256)
c = CentralCharge(β = big"1.2" + big"0.1" * im)
c_s = CentralCharge(β = c.β / sqrt(big(2)))

V = Field(c, P = big"0.53" + big"0.11"*im, diagonal = true)
V2 = Field(c, P = big"0.43" + big"0.31"*im, diagonal = true)
V1 = Field(c, P = big"0.71" + big"1.03"*im, diagonal = true)
d = (V1,)

Nmax = 40
co = Correlation(V1, Nmax)

τ = big"0.3" + big"2" * im # τ in H/PSL_2(ZZ)
τ_cache = PositionCache(τ, co, :s)
q = qfromτ(τ)
q_s = exp(im * (π * τ))
x = xfromq(exp(im * (π * τ)))

V_s = Field(c_s, P = sqrt(big(2)) * V.dims[:left].P)
V2_s = Field(c_s, P = sqrt(big(2)) * V2.dims[:left].P)
V1_s = Field(c_s, P = V1.dims[:left].P / sqrt(big(2)))
V_kac = Field(c_s, r = 0, s = 1 // 2)
Vs = (V_kac, V1_s, V_kac, V_kac)

co_s = Correlation(Vs..., Nmax)

@testset "Chiral" begin
    # Match the first residue
    co_l = co[:left]
    co_l_s = co_s[:left]


    co_l_s._Rmn[:s][2, 1] * 16^2 / 2
    co_l._Rmn[:s][1, 1]

    @test isapprox(co_l_s._Rmn[:s][2, 1] * 16^2 / 2, co_l._Rmn[:s][1, 1])

    @test isapprox(co_l_s._CNmn[:s][2, 2, 1] * 16^2 / 2, co_l._CNmn[:s][1, 1, 1])

    @test isapprox(co_l_s._Rmn[:s][4, 1] * 16^4 / 2, co_l._Rmn[:s][2, 1])

    @test isapprox(co_l_s._Rmn[:s][2, 2] * 16^4 / 2, co_l._Rmn[:s][1, 2])

    @test isapprox(co_l_s._CNmn[:s][6, 2, 2] * 16^6 / 2, co_l._CNmn[:s][3, 1, 2])

    q = τ_cache.q
    @test qfromx(xfromq(q)) ≈ q
    @test qfromx(x)^2 ≈ q

    b = Block(co, :s, V, :left)
    b_s = Block(co_s, :s, V_s, :left)

    b2 = Block(co, :s, V2, :left)
    b2_s = Block(co_s, :s, V2_s, :left)

    h = evaluate_series(b, τ_cache)
    @test h ≈ Base.evalpoly(τ_cache.q, b._coefficients)

    # q_s = 16BootstrapVirasoro.qfromx(x)

    # @test q_s ≈ 16 * q

    h_s = evaluate_series(b_s, x)
    @test isapprox(h, h_s, rtol = 1e-14)

    F = evaluate(b, τ)
    F_s = evaluate(b_s, x)

    F2 = evaluate(b2, τ)
    F2_s = evaluate(b2_s, x)

    @test isapprox(
        F / BootstrapVirasoro.total_prefactor(b, τ_cache),
        F_s / BootstrapVirasoro.total_prefactor(b_s, x),
        rtol = 1e-16,
    )

    @test isapprox(
        F / F2 * 16^(2 * (V.dims[:left].P^2 - V2.dims[:left].P^2)),
        F_s / F2_s,
        rtol = 1e-15,
    ) # up to P-indep prefactors 16^(2P^2)*F_P^t = F_{\sqrt{2} P}^s
end

@testset "Fder" begin
    setprecision(BigFloat, 256)
    ϵ = 1e-25

    V0 = Field(c, :P, big"0.5", diagonal = true)
    Vp = Field(c, :P, big"0.5" + ϵ, diagonal = true)
    Vm = Field(c, :P, big"0.5" - ϵ, diagonal = true)

    b = Block(co, :s, V0, :left, 40, der = true)
    b_noder = Block(co, :s, V0, :left, 40)
    b_p = Block(co, :s, Vp, :left, 40)
    b_m = Block(co, :s, Vm, :left, 40)

    @test isapprox(
        evaluate(b, 0.3 + 0.4im),
        evaluate(b_noder, 0.3 + 0.4im),
        rtol = 1e-15,
    )

    function eval_series_der(τ)
        series_der = evaluate_series_der(b, τ)
        byhand = (evaluate_series(b_p, τ) - evaluate_series(b_m, τ)) / (2 * ϵ)
        series_der - byhand
    end

    @test abs(eval_series_der(big"0.3" + big"0.4" * im)) < big"1e-30"
    @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-30"

    function eval_block_der(τ)
        block_der = evaluate_der(b, τ)
        byhand = (evaluate(b_p, τ) - evaluate(b_m, τ)) / (2 * ϵ)
        block_der - byhand
    end

    import BootstrapVirasoro: etaDedekind

    byhand = (evaluate(b_p, τ) - evaluate(b_m, τ)) / (2 * ϵ)
    h = evaluate_series(b, τ)
    hprime_byhand = (evaluate_series(b_p, τ) - evaluate_series(b_m, τ)) / (2 * ϵ)

    @test isapprox(hprime_byhand, evaluate_series_der(b, τ))

    @test isapprox(
        evaluate(b, τ) / evaluate_series(b, τ),
        q^(V0.dims[:left].δ) / etaDedekind(τ),
        rtol = 1e-50,
    )

    @test isapprox(
        evaluate_der(b, τ),
        q^(V0.dims[:left].δ) / etaDedekind(τ) *
        (hprime_byhand + 2 * V0.dims[:left].P * log(q) * h),
        rtol = 1e-30,
    )

    test_vals =
        [big"0.3" + big"0.4" * im, big"0.6" + big"1.2" * im, big"10" + big"0.01" * im]
    for val in test_vals
        @test abs(eval_block_der(val)) < big"1e-23"
    end
end

@testset "Freg" begin
    V0 = Field(c, r = 2, s = 1)
    ϵ = 1e-40
    δ0ϵ = V0.dims[:left].δ + ϵ
    V0ϵ = Field(c, diagonal = true, δ = δ0ϵ)
    V0_s = Field(c_s, r = 2 * V0.r, s = V0.s)
    δ0_sϵ = V0_s.dims[:left].δ + ϵ
    V0_sϵ = Field(c_s, diagonal = true, δ = δ0_sϵ)
    b = Block(co, :s, V0, Nmax)
    bϵ = Block(co, :s, V0ϵ, Nmax)
    b_s = Block(co_s, :s, V0_s, Nmax)
    b_sϵ = Block(co_s, :s, V0ϵ, Nmax)

    lr_cache = BootstrapVirasoro.LRPositionCache(τ, b)
    xlr_cache = BootstrapVirasoro.LRPositionCache(x, b_s)
    @test isapprox(
        evaluate(bϵ[:left], τ),
        co._Rmn[:left][:s][V0.r, V0.s] / ϵ *
        BootstrapVirasoro.evaluate_lr_op(b, lr_cache)[:left] + evaluate(b[:left], τ),
        rtol = 1e-40,
    ) # F_{Prs + ϵ} = R/ϵ F_{Pr,-s} + F^reg

    import BootstrapVirasoro: ell, etaDedekind, conj_q
    import BootstrapVirasoro.total_prefactor as prefac
    import BootstrapVirasoro: evaluate_series

    @test isapprox(
        evaluate_series(b[:left], τ_cache),
        evaluate_series(b_s[:left], x),
        rtol = 1e-40,
    ) # regularized series

    Freg = evaluate(b[:left], τ)
    Rrs = co._Rmn[:left][:s][V0.r, V0.s]
    Fminus = BootstrapVirasoro.evaluate_lr_op(b, lr_cache)[:left]
    Freg_s = evaluate(b_s[:left], x)
    Rrs_s = co_s._Rmn[:left][:s][V0_s.r, V0_s.s]
    Fminus_s = BootstrapVirasoro.evaluate_lr_op(b_s, xlr_cache)[:left]

    @test isapprox(
        Freg_s / prefac(b_s[:left], x),
        (Freg + 8log(big"2") * Rrs * Fminus) / prefac(b[:left], τ),
        rtol = 1e-40,
    ) # regularized block

    @test isapprox(
        (Freg_s - 4log(big"2") * Rrs_s * Fminus_s) / prefac(b_s[:left], x),
        Freg / prefac(b[:left], τ),
        rtol = 1e-40,
    ) # regularized block
end

@testset "Flog" begin
    import BootstrapVirasoro: ell, etaDedekind, conj_q, total_prefactor

    V = Field(c, r = 2, s = 1)
    V_s = Field(c_s, r = 2 * V.r, s = V.s)

    b = Block(co, :s, V, Nmax)
    b_op = Block(co, :s, Field(V.c, r = V.r, s = -V.s), Nmax)
    b_s = Block(co_s, :s, V_s, Nmax)
    b_op_s = Block(co_s, :s, Field(V_s.c, r = V_s.r, s = -V_s.s), Nmax)

    P23 = V.dims[:left].P
    P43_s = V_s.dims[:left].P

    @test isapprox(
        evaluate_series(b[:left], τ_cache),
        evaluate_series(b_s[:left], x),
        rtol = 1e-40,
    ) # regularized series match

    missing_terms = b[:left]._missing_terms
    missing_term = log(q) * Base.evalpoly(q, missing_terms)

    @test isapprox(
        evaluate(b[:left], τ) / total_prefactor(b[:left], τ),
        evaluate_series(b[:left], τ_cache) + missing_term,
        rtol = 1e-40,
    ) # the regularised block is as expected

    missing_terms_s = b_s[:left]._missing_terms
    missing_term_s = log(16q_s) * Base.evalpoly(16q_s, missing_terms_s)

    @test 2 * Base.evalpoly(q, missing_terms) ≈ Base.evalpoly(16q_s, missing_terms_s)

    @test isapprox(log(q_s) / log(16q_s) * missing_term_s, missing_term, rtol = 1e-40)

    @test isapprox(
        evaluate(b_s[:left], x) / total_prefactor(b_s[:left], x),
        evaluate_series(b_s[:left], x) + missing_term_s,
        rtol = 1e-40,
    ) # the regularised block is as expected

    @test isapprox(
        evaluate(b[:left], τ) / total_prefactor(b[:left], τ) - missing_term,
        evaluate(b_s[:left], x) / total_prefactor(b_s[:left], x) - missing_term_s,
        rtol = 1e-40,
    )
end

import BootstrapVirasoro: ell, etaDedekind, conj_q, total_prefactor

V = Field(c, r = 2, s = 1)
V2 = Field(c, r = 4, s = 4)
V_s = Field(c_s, r = 2 * V.r, s = V.s)
V2_s = Field(c_s, r = 2 * V2.r, s = V2.s)

b = Block(co, :s, V, Nmax)
b2 = Block(co, :s, V2, Nmax)
b_s = Block(co_s, :s, V_s, Nmax)
b2_s = Block(co_s, :s, V2_s, Nmax)

pref = total_prefactor(b[:left], τ) * total_prefactor(b[:right], conj_q(τ, b.corr))
pref_s =
    total_prefactor(b_s[:left], x) * total_prefactor(b_s[:right], conj_q(x, b_s.corr))
pref2 =
    total_prefactor(b2[:left], τ) * total_prefactor(b2[:right], conj_q(τ, b.corr))
pref2_s =
    total_prefactor(b2_s[:left], x) *
    total_prefactor(b2_s[:right], conj_q(x, b_s.corr))

@testset "ell" begin
    @test isapprox(
        ell(b) / sqrt(big"2") + 16 * log(big"2") * inv(c_s.β) * V.s,
        ell(b_s),
    ) # ell/sqrt(2) + 16 log(2) \beta'^{-1} s = ell^{S^2}

    @test isapprox(
        ell(b2) / sqrt(big"2") + 16 * log(big"2") * inv(c_s.β) * V2.s,
        ell(b2_s),
    ) # ell/sqrt(2) + 16 log(2) \beta'^{-1} s = ell^{S^2}
end

@testset "Full log block" begin
    @test isapprox(evaluate(b, τ) / pref, evaluate(b_s, x) / pref_s)

    @test isapprox(evaluate(b2, τ) / pref2, evaluate(b2_s, x) / pref2_s)
end

@testset "Interchiral" begin
    import BootstrapVirasoro: shift_D, conj_q, total_prefactor

    P = big"0.53" + big"0.11" * im
    V = Field(c, P = P, diagonal = true)
    V_s = Field(c_s, P = sqrt(big"2") * P, diagonal = true)

    s = shift_D(co.fields, V)
    s_s = shift_D(co_s.fields, V_s)

    s_s / s
    16^(8 / c.β * (V.dims[:left].P - 1/2/c.β))

    @test isapprox(s_s / s / 16^(8 / c.β * (V.dims[:left].P - 1 / 2 / c.β)), 1) # shift(D^S2) = shift(D) * 16^(-8β^-1 (P - β^-1/2))

    b = Block(co, :s, V, interchiral = true, Δmax = 10)
    prefactor =
        BootstrapVirasoro.PositionCache(2τ, b.blocks[1][:left]).prefactor *
        BootstrapVirasoro.PositionCache(
            conj_q(2τ, b.blocks[1].corr),
            b.blocks[1][:right],
        ).prefactor

    b_s = Block(co_s, :s, V_s, interchiral = true, Δmax = 10)
    prefactor_s =
        BootstrapVirasoro.PositionCache(x, b_s.blocks[1][:left]).prefactor *
        BootstrapVirasoro.PositionCache(
            conj_q(x, b_s.blocks[1].corr),
            b_s.blocks[1][:right],
        ).prefactor

    for i = 1:3
        @test isapprox(
            evaluate(b.blocks[i], 2τ) / prefactor,
            evaluate(b_s.blocks[i], x) / prefactor_s *
            16^(-4b.fields[i].dims[:left].P^2),
            atol = 1e-20,
        )
    end

    @test isapprox(
        (evaluate(b.blocks[1], 2τ) + evaluate(b.blocks[2], 2τ) * s) / prefactor,
        (evaluate(b_s.blocks[1], x) + evaluate(b_s.blocks[2], x) * s_s) /
        prefactor_s * 16^(-4P^2),
        rtol = 1e-20,
    )

    @test isapprox(
        evaluate(b, 2τ) / prefactor,
        evaluate(b_s, x) / prefactor_s * 16^(-4P^2),
        rtol = 1e-20,
    )
end
