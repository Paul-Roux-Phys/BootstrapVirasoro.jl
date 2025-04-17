using BootstrapVirasoro: series_argument, blockprefactor_chiral, islogarithmic

setprecision(BigFloat, 256)
c = CentralCharge(β=big"1.2" + big"0.1" * 1im)
c_s = CentralCharge(β=c.β / sqrt(big(2)))

V = Field(c, P=0.53 + 0.11im, diagonal=true)
V2 = Field(c, P=0.43 + 0.31im, diagonal=true)
V1 = Field(c, P=0.71 + 1.03im, diagonal=true)
d = (V1,)

Nmax = 40
co = Correlation(V1, Nmax)

τ = big"0.3" + big"2" * im # τ in H/PSL_2(ZZ)
q = exp(im * big(π) * τ)
x = xfromq(q)

V_s = Field(c_s, P=sqrt(2) * V.P[:left])
V2_s = Field(c_s, P=sqrt(2) * V2.P[:left])
V1_s = Field(c_s, P=V1.P[:left] / sqrt(big(2)))
V_kac = Field(c_s, r=0, s=1 // 2)
Vs = (V_kac, V1_s, V_kac, V_kac)

co_s = Correlation(Vs..., Nmax)

@testset "Chiral" begin
    # Match the first residue
    co_l = co[:left]
    co_l_s = co_s[:left]
    @test isapprox(
        co_l_s._Rmn[:s][(2, 1)] * 16^2 / 2,
        co_l._Rmn[:τ][(1, 1)]
    )

    @test isapprox(
        co_l_s._CNmn[:s][(2, 2, 1)] * 16^2 / 2,
        co_l._CNmn[:τ][(1, 1, 1)]
    )

    @test isapprox(
        co_l_s._Rmn[:s][(4, 1)] * 16^4 / 2,
        co_l._Rmn[:τ][(2, 1)]
    )

    @test isapprox(
        co_l_s._Rmn[:s][(2, 2)] * 16^4 / 2,
        co_l._Rmn[:τ][(1, 2)]
    )

    @test isapprox(
        co_l_s._CNmn[:s][(6, 2, 2)] * 16^6 / 2,
        co_l._CNmn[:τ][(3, 1, 2)]
    )

    @test isapprox(qfromx(xfromq(q)), q)
    @test isapprox(qfromx(x), q)

    b = Block(co, :τ, V, :left)
    b_s = Block(co_s, :s, V_s, :left)

    b2 = Block(co, :τ, V2, :left)
    b2_s = Block(co_s, :s, V2_s, :left)

    y_s = BootstrapVirasoro.get_position(x, Vs, b_s)
    q_s = BootstrapVirasoro.series_argument(y_s, b_s)

    @test isapprox(q_s, 16 * q)

    h = evaluate_series(b, q^2)
    h_s = evaluate_series(b_s, q_s)
    @test isapprox(h, h_s, rtol=1e-14)

    F = evaluate(b, 2τ)
    F_s = evaluate(b_s, x)

    F2 = evaluate(b2, 2τ)
    F2_s = evaluate(b2_s, x)

    @test isapprox(
        F / BootstrapVirasoro.blockprefactor_chiral(b, 2τ),
        F_s / BootstrapVirasoro.blockprefactor_chiral(b_s, x),
        rtol=1e-16
    )

    @test isapprox(
        F / F2 * 16^(2 * (V.P[:left]^2 - V2.P[:left]^2)), F_s / F2_s, rtol=1e-15
    ) # up to P-indep prefactors 16^(2P^2)*F_P^t = F_{\sqrt{2} P}^s
end

@testset "Non Chiral" begin
    @testset "Fder" begin
        setprecision(BigFloat, 256)
        ϵ = 1e-25

        V0 = Field(c, :P, big"0.5", diagonal=true)
        Vp = Field(c, :P, big"0.5" + ϵ, diagonal=true)
        Vm = Field(c, :P, big"0.5" - ϵ, diagonal=true)

        b = Block(co, :τ, V0, :left, 40, der=true)
        b_noder = Block(co, :τ, V0, :left, 40)
        b_p = Block(co, :τ, Vp, :left, 40)
        b_m = Block(co, :τ, Vm, :left, 40)

        @test isapprox(
            evaluate(b, 0.3 + 0.4im),
            evaluate(b_noder, 0.3 + 0.4im),
            rtol=1e-15
        )

        function eval_series_der(τ)
            q = exp(im * π * τ)
            series_der = evaluate_series(b, q, der=true)
            byhand = (evaluate_series(b_p, q) - evaluate_series(b_m, q)) / (2 * ϵ)
            series_der - byhand
        end

        @test abs(eval_series_der(big"0.3" + big"0.4" * im)) < big"1e-45"
        @test abs(eval_series_der(big"0.01")) < big"1e-30"
        @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-30"

        function eval_block_der(τ)
            block_der = evaluate(b, τ, der=true)
            byhand = (evaluate(b_p, τ) - evaluate(b_m, τ)) / (2 * ϵ)
            block_der - byhand
        end

        import BootstrapVirasoro: etaDedekind

        byhand = (evaluate(b_p, τ) - evaluate(b_m, τ)) / (2 * ϵ)
        h = evaluate_series(b, q)
        hprime_byhand = (evaluate_series(b_p, q) - evaluate_series(b_m, q)) / (2 * ϵ)

        @test isapprox(
            hprime_byhand,
            evaluate_series(b, q, der=true)
        )

        @test isapprox(
            evaluate(b, τ) / evaluate_series(b, q),
            q^(V0.δ[:left]) / etaDedekind(τ),
            rtol=1e-50
        )

        @test isapprox(
            evaluate(b, τ, der=true),
            q^(V0.δ[:left]) / etaDedekind(τ) * (hprime_byhand + 2 * V0.P[:left] * log(q) * h),
            rtol=1e-30
        )

        test_vals = [big"0.3" + big"0.4" * im, big"0.6" + big"1.2" * im, big"10" + big"0.01" * im]
        for val in test_vals
            @test abs(eval_block_der(val)) < big"1e-23"
        end
    end

    @testset "Freg" begin
        V0 = Field(c, r=2, s=1)
        ϵ = 1e-40
        δ0ϵ = V0.δ[:left] + ϵ
        V0ϵ = Field(c, diagonal=true, δ=δ0ϵ)
        V0_s = Field(c_s, r=2 * V0.r, s=V0.s)
        δ0_sϵ = V0_s.δ[:left] + ϵ
        V0_sϵ = Field(c_s, diagonal=true, δ=δ0_sϵ)
        b = Block(co, :τ, V0, Nmax)
        bϵ = Block(co, :τ, V0ϵ, Nmax)
        b_s = Block(co_s, :s, V0_s, Nmax)
        b_sϵ = Block(co_s, :s, V0ϵ, Nmax)

        @test isapprox(
            evaluate(bϵ, 2τ, :left),
            co._Rmn[:left][:τ][(V0.r, V0.s)] / ϵ * evaluate(b, 2τ, :left, op=true) + evaluate(b, 2τ, :left),
            rtol=1e-40
        ) # F_{Prs + ϵ} = R/ϵ F_{Pr,-s} + F^reg

        import BootstrapVirasoro: ell, etaDedekind, conj_q
        import BootstrapVirasoro.blockprefactor_chiral as prefac
        import BootstrapVirasoro: evaluate_series

        @test isapprox(
            evaluate_series(b[:left], q^2),
            evaluate_series(b_s[:left], 16q),
            rtol=1e-40
        ) # regularized series

        Freg = evaluate(b, 2τ, :left)
        Rrs = co._Rmn[:left][:τ][(V0.r, V0.s)]
        Fminus = evaluate(b, 2τ, :left, op=true)
        Freg_s = evaluate(b_s, x, :left)
        Rrs_s = co_s._Rmn[:left][:s][(V0_s.r, V0_s.s)]
        Fminus_s = evaluate(b_s, x, :left, op=true)

        @test isapprox(
            Freg_s / prefac(b_s[:left], x),
            (Freg + 8log(big"2") * Rrs * Fminus) / prefac(b[:left], 2τ),
            rtol=1e-40
        ) # regularized block

        @test isapprox(
            (Freg_s - 4log(big"2") * Rrs_s * Fminus_s) / prefac(b_s[:left], x),
            Freg / prefac(b[:left], 2τ),
            rtol=1e-40
        ) # regularized block

        @test isapprox(
            evaluate(b, 2τ, :right) / prefac(b[:right], conj_q(2τ, b)),
            evaluate(b_s, x, :right) / prefac(b_s[:right], conj_q(x, b_s)),
            rtol=1e-40
        ) # non regularized block
    end

    @testset "Flog" begin
        import BootstrapVirasoro: ell, etaDedekind, conj_q
        import BootstrapVirasoro.blockprefactor_chiral as prefac

        V = Field(c, r=2, s=1)
        V_s = Field(c_s, r=2 * V.r, s=V.s)

        b = Block(co, :τ, V, Nmax)
        b_op = Block(co, :τ, Field(V.c, r=V.r, s=-V.s), Nmax)
        b_s = Block(co_s, :s, V_s, Nmax)
        b_op_s = Block(co_s, :s, Field(V_s.c, r=V_s.r, s=-V_s.s), Nmax)

        P23 = V.P[:left]
        P43_s = V_s.P[:left]

        @test isapprox(
            evaluate_series(b[:left], q^2),
            evaluate_series(b_s[:left], 16q),
            rtol=1e-40
        ) # regularized series match

        missing_terms = [
            (N, V.r, V.s) in keys(b[:left]._CNmn) ? b[:left]._CNmn[(N, V.r, V.s)] : zero(x)
            for N in 0:Nmax
        ]
        missing_term = log(q^2) * evalpoly(q^2, missing_terms)

        @test isapprox(
            evaluate(b, 2τ, :left) / prefac(b[:left], 2τ),
            evaluate_series(b[:left], q^2) + missing_term,
            rtol=1e-40
        ) # the regularised block is as expected

        missing_terms_s = [
            (N, V_s.r, V_s.s) in keys(b_s[:left]._CNmn) ? b_s[:left]._CNmn[(N, V_s.r, V_s.s)] : zero(x)
            for N in 0:Nmax
        ]
        missing_term_s = log(16q) * evalpoly(16q, missing_terms_s)

        @test isapprox(
            log(q) / log(16q) * missing_term_s, missing_term, rtol=1e-40
        )

        @test isapprox(
            evaluate(b_s, x, :left) / prefac(b_s[:left], x),
            evaluate_series(b_s[:left], 16q) + missing_term_s,
            rtol=1e-40
        ) # the regularised block is as expected

        @test isapprox(
            evaluate(b, 2τ, :left) / prefac(b[:left], 2τ) - missing_term,
            evaluate(b_s, x, :left) / prefac(b_s[:left], x) - missing_term_s,
            rtol=1e-40
        )

        logterm1 = evaluate(b, 2τ, debug=true)[1]
        pref = prefac(b[:left], 2τ) * prefac(b[:right], conj_q(2τ, b))
        logterm1_s = evaluate(b_s, x, debug=true)[1]
        R2rs_s = co_s._Rmn[:left][:s][(V_s.r, V_s.s)]
        pref_s = prefac(b_s[:left], x) * prefac(b_s[:right], conj_q(x, b_s))
        F2rms_s = evaluate(b_s, x, :right)
        F2Rmsbar_s = evaluate(b_s, conj_q(x, b_s), :right)
        modF2sq = F2rms_s * F2Rmsbar_s
        Prs = V.P[:left]

        @test isapprox(
            logterm1 / pref,
            (logterm1_s + 4log(big"2") * R2rs_s / c.β * V_s.s / Prs * modF2sq) / pref_s
        ) # the terms F^log obey the correct relation
    end

    import BootstrapVirasoro: ell, etaDedekind, conj_q
    import BootstrapVirasoro.blockprefactor_chiral as prefac

    V = Field(c, r=2, s=1)
    V2 = Field(c, r=4, s=4)
    V_s = Field(c_s, r=2 * V.r, s=V.s)
    V2_s = Field(c_s, r=2 * V2.r, s=V2.s)

    b = Block(co, :τ, V, Nmax)
    b2 = Block(co, :τ, V2, Nmax)
    b_s = Block(co_s, :s, V_s, Nmax)
    b2_s = Block(co_s, :s, V2_s, Nmax)

    pref = prefac(b[:left], 2τ) * prefac(b[:right], conj_q(2τ, b))
    pref_s = prefac(b_s[:left], x) * prefac(b_s[:right], conj_q(x, b_s))
    pref2 = prefac(b2[:left], 2τ) * prefac(b2[:right], conj_q(2τ, b))
    pref2_s = prefac(b2_s[:left], x) * prefac(b2_s[:right], conj_q(x, b_s))

    setprecision(BigFloat, 256)
    println(ell(b2) / sqrt(big"2") + 16 * log(big"2") * inv(c_s.β) * V2.s)
    println(ell(b2_s))

    @testset "ell" begin
        @test isapprox(
            ell(b) / sqrt(big"2") + 16 * log(big"2") * inv(c_s.β) * V.s, ell(b_s)
        ) # ell/sqrt(2) + 16 log(2) \beta'^{-1} s = ell^{S^2}

        @test isapprox(
            ell(b2)/sqrt(big"2") + 16*log(big"2") * inv(c_s.β) * V2.s, ell(b2_s)
        ) # ell/sqrt(2) + 16 log(2) \beta'^{-1} s = ell^{S^2}
    end

    @testset "Full log block" begin
        @test isapprox(
            evaluate(b, 2τ) / pref, evaluate(b_s, x) / pref_s
        )

        @test isapprox(
            evaluate(b2, 2τ) / pref2, evaluate(b2_s, x) / pref2_s
        )
    end
end
