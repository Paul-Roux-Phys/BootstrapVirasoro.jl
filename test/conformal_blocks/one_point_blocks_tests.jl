using BootstrapVirasoro: series_argument

c = CentralCharge(β=big"1.2" + big"0.1" * 1im)
c_s = CentralCharge(β=c.β / sqrt(big(2)))

V = Field(c, P=0.53 + 0.11im, diagonal=true)
V1 = Field(c, P=0.71 + 1.03im, diagonal=true)
d = (V1,)

Nmax = 26
co = Correlation(V1, Nmax)

@testset "Relation to sphere 4pt" begin
    V_s = Field(c_s, P=sqrt(2) * V.P[:left])
    V1_s = Field(c_s, P=V1.P[:left] / sqrt(big(2)))
    V_kac = Field(c_s, r=0, s=1 // 2)
    Vs = (V_kac, V1_s, V_kac, V_kac)

    co_s = Correlation(Vs..., Nmax)

    @testset "Chiral" begin
        # Match the first residue
        co_l = co[:left]
        co_l_s = co_s[:left]
        @test isapprox(
            co_l_s._Rmn[:s][(2, 1)] * 2^7,
            co_l._Rmn[:τ][(1, 1)]
        )

        @test isapprox(
            co_l_s._CNmn[:s][(2, 2, 1)] * 2^7,
            co_l._CNmn[:τ][(1, 1, 1)]
        )

        @test isapprox(
            co_l_s._Rmn[:s][(4, 1)] * 2^15,
            co_l._Rmn[:τ][(2, 1)]
        )

        @test isapprox(
            co_l_s._Rmn[:s][(2, 2)] * 2^15,
            co_l._Rmn[:τ][(1, 2)]
        )

        τ = 0.3 + 2im # τ in H/PSL_2(ZZ)
        q = exp(im * big(π) * τ)
        x = xfromq(q)

        @test isapprox(qfromx(xfromq(q)), q)
        @test isapprox(qfromx(x), q)

        y_s = BootstrapVirasoro.get_position(x, Vs, b_s)
        q_s = BootstrapVirasoro.series_argument(y_s, b_s)
        @test isapprox(q_s, 16 * q)

        b = Block(co, :τ, V, :left)
        b_s = Block(co_s, :s, V_s, :left)

        @test isapprox(
            evaluate_series(b, q^2),
            evaluate_series(b_s, q_s),
            rtol=1e-14
        )

        F = evaluate(b, 2 * τ)
        F_s = evaluate(b_s, x)

        @test isapprox(
            F / BootstrapVirasoro.blockprefactor_chiral(Tuple(D.dims[:left] for D in d), b, 2 * τ), # F^T(q^2)/prefactor
            F_s / BootstrapVirasoro.blockprefactor_chiral(Tuple(v.dims[:left] for v in Vs), b_s, x),  # F^S(q)/prefactor
            rtol=1e-16
        )
    end

    @testset "Non Chiral" begin
        V = Field(c, r=2, s=3)
        x = big"0.3" + big"0.1" * im

        @testset "Derivatives" begin
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
            @test abs(eval_series_der(big"0.01")) < big"1e-35"
            @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-35"

            function eval_block_der(τ)
                block_der = evaluate(b, τ, der=true)
                byhand = (evaluate(b_p, τ) - evaluate(b_m, τ)) / (2 * ϵ)
                block_der - byhand
            end

            import BootstrapVirasoro: etaDedekind

            τ = big"0.3" + big"0.5" * im
            q = exp(im * (π * τ))
            byhand = (evaluate(b_p, τ) - evaluate(b_m, τ)) / (2 * ϵ)
            h = evaluate_series(b, q)
            hprime_byhand = (evaluate_series(b_p, q) - evaluate_series(b_m, q)) / (2 * ϵ)

            @test isapprox(
                hprime_byhand,
                evaluate_series(b, q, der=true)
            )

            @test isapprox(
                evaluate(b_p, τ) / evaluate_series(b_p, q),
                q^(V0.δ[:left]) / etaDedekind(τ),
                rtol=1e-14
            )

            evaluate(b_p, τ)
            evaluate(b_m, τ)

            @test isapprox(
                evaluate(b, τ, der=true),
                q^(V0.δ[:left]) / etaDedekind(τ) * (hprime_byhand + 2 * V0.P[:left] * log(q) * h),
                rtol=1e-40
            )

            @test abs(eval_block_der(big"0.3" + big"0.4" * im)) < big"1e-25"
            @test abs(eval_block_der(big"0.6" + big"1.2" * im)) < big"1e-25"
            @test abs(eval_block_der(big"10" + big"0.01" * im)) < big"1e-25"
        end

        @testset "Logarithmic blocks" begin
            setprecision(BigFloat, 256)
            ϵ = 1e-25
            V0 = Field(c, r=2, s=3)
            b = Block(co, :τ, V0, 40)

            import BootstrapVirasoro: ell
            l = ell(b, 2, 1)

            F = evaluate(b, 2 * τ)
            F_s = evaluate(b_s, x)

            @test true
        end
    end
end
