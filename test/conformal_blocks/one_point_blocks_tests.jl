using BootstrapVirasoro: series_argument
@testset "Relation to sphere four-point blocks" begin
    c_t = CentralCharge(β=1.2 + 0.1 * 1im)
    c_s = CentralCharge(β=c_t.β / sqrt(2))

    V_t = ConformalDimension(c_t, P=0.53 + 0.11im)
    V1_t = ConformalDimension(c_t, P=0.71 + 1.03im)
    d = (V1_t,)

    V_s = ConformalDimension(c_s, P=sqrt(2) * V_t.P)
    V1_s = ConformalDimension(c_s, P=V1_t.P / sqrt(2))
    V_kac = ConformalDimension(c_s, r=0, s=1 // 2)
    ds = (V_kac, V1_s, V_kac, V_kac)

    co_t = Correlation(V1_t, 1)
    b_t = Block(co_t, :τ, V_t)

    co_s = Correlation(ds..., 2)
    b_s = Block(co_s, :s, V_s)

    # Match the first residue
    @test isapprox(
        co_s._Rmn[:s][(2, 1)] * 2^7,
        co_t._Rmn[:τ][(1, 1)]
    )

    @test isapprox(
        co_s._CNmn[:s][(2, 2, 1)] * 2^7,
        co_t._CNmn[:τ][(1, 1, 1)]
    )

    τ = 0.3 + 2im # τ in H/PSL_2(Z)
    q = exp(im * big(π) * τ)
    x = xfromq(q)

    @test isapprox(qfromx(xfromq(q)), q)
    @test isapprox(qfromx(x), q)

    y_s = BootstrapVirasoro.get_position(x, ds, b_s)
    q_s = BootstrapVirasoro.series_argument(y_s, b_s)
    @test isapprox(q_s, 16 * q)

    @test isapprox(
        evaluate_series(b_t, q^2),
        evaluate_series(b_s, q_s),
        rtol = 1e-12
    )

    evaluate_series(b_t, q^2)

    F_t = evaluate(b_t, τ)

    F_t/BootstrapVirasoro.blockprefactor_chiral(d, b_t, τ)
    
    F_s = evaluate(b_s, x)

    @test isapprox(
        F_t / BootstrapVirasoro.blockprefactor_chiral(d, b_t, τ),
        F_s / BootstrapVirasoro.blockprefactor_chiral(ds, b_s, x),
        atol=1e-15
    )
end
