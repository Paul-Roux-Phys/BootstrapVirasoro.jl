import BootstrapVirasoro:
    LRPosCache,
    PosCache,
    eval_series,
    etaDedekind,
    ell,
    conj_q,
    xfromq,
    qfromx,
    qfromτ,
    τfromx,
    total_prefactor,
    shift_D
import BootstrapVirasoro.total_prefactor as prefac

setprecision(BigFloat, 256)
c = CentralCharge(β = big"1.2" + big"0.1" * im)
c_S2 = CentralCharge(β = c.β / sqrt(big(2)))

V = Field(c, P = big"0.53" + big"0.11" * im, diagonal = true)
V1 = Field(c, P = big"0.71" + big"1.03" * im, diagonal = true)
V2 = Field(c, P = big"0.43" + big"0.31" * im, diagonal = true)

Δmax = 40
co = Correlation(V1, Δmax)

τ = big"0.3" + big"2" * im # τ in H/PSL_2(ZZ)
τ_cache = LRPosCache(τ, co, :s)
q = qfromτ(τ)
q_S2 = exp(im * (π * τ))
x = xfromq(exp(im * (π * τ)))

V_S2 = Field(c_S2, P = sqrt(big(2)) * V.dims.left.P)
V2_S2 = Field(c_S2, P = sqrt(big(2)) * V2.dims.left.P)
V1_S2 = Field(c_S2, P = V1.dims.left.P / sqrt(big(2)))
V_kac = Field(c_S2, r = 0, s = 1 // 2)
ind_S2 = [(0, 1 // 2), (0, V1.s / 2), (0, 1 // 2), (0, 1 // 2)]

co_S2 = Correlation([Field(c_S2, r = r, s = s) for (r, s) in ind_S2], Δmax)

@testset "Chiral" begin
    # Match the first residue
    co_l = co[:left]
    co_l_S2 = co_S2[:left]

    @test co_l_S2.Rmn[:s][2, 1] * 16^2 / 2 ≈ co_l.Rmn[1, 1]
    @test co_l_S2.CNmn[:s][2, 2, 1] * 16^2 / 2 ≈ co_l.CNmn[1, 1, 1]
    @test co_l_S2.Rmn[:s][4, 1] * 16^4 / 2 ≈ co_l.Rmn[2, 1]
    @test co_l_S2.Rmn[:s][2, 2] * 16^4 / 2 ≈ co_l.Rmn[1, 2]
    @test co_l_S2.CNmn[:s][6, 2, 2] * 16^6 / 2 ≈ co_l.CNmn[3, 1, 2]

    q = τ_cache.left.q
    @test qfromx(xfromq(q)) ≈ q
    @test qfromx(x)^2 ≈ q

    b = Block(co_l, V.dims.left)
    b_S2 = Block(co_l_S2, V_S2.dims.left)

    b2 = Block(co_l,V2.dims.left)
    b2_S2 = Block(co_l_S2, V2_S2.dims.left)

    h = eval_series(b, τ_cache.left)
    @test h ≈ Base.evalpoly(q, b.coeffs)

    h_S2 = eval_series(b_S2, x)
    @test h ≈ h_S2

    F = b(τ)
    F_S2 = b_S2(x)

    F2 = b2(τ)
    F2_S2 = b2_S2(x)

    @test F / BootstrapVirasoro.total_prefactor(b, τ_cache.left) ≈
          F_S2 / BootstrapVirasoro.total_prefactor(b_S2, x)

    @test F / F2 * 16^(2 * (V.dims[:left].P^2 - V2.dims[:left].P^2)) ≈ F_S2 / F2_S2
    # up to P-indep prefactors 16^(2P^2)*F_P^t = F_{\sqrt{2} P}^s
end

@testset "Fder" begin
    ϵ = 1e-25

    V0 = Field(c, P = big"0.5")
    Vp = Field(c, P = big"0.5" + ϵ)
    Vm = Field(c, P = big"0.5" - ϵ)

    co_l = co[:left]
    b = Block(co_l, :s, V0.dims.left, 40, der = true)
    b_noder = Block(co_l, :s, V0.dims.left, 40)
    b_p = Block(co_l, :s, Vp.dims.left, 40)
    b_m = Block(co_l, :s, Vm.dims.left, 40)

    @test isapprox(b(0.3 + 0.4im), b_noder(0.3 + 0.4im), rtol = 1e-15)

    function eval_series_der(τ)
        series_der = BootstrapVirasoro.eval_series_der(b, τ)
        byhand =
            (
                BootstrapVirasoro.eval_series(b_p, τ) -
                BootstrapVirasoro.eval_series(b_m, τ)
            ) / (2 * ϵ)
        series_der - byhand
    end

    @test abs(eval_series_der(big"0.3" + big"0.4" * im)) < big"1e-30"
    @test abs(eval_series_der(big"10" + big"0.01" * im)) < big"1e-30"

    function eval_block_der(τ)
        block_der = BootstrapVirasoro.eval_der(b, τ)
        byhand =
            (BootstrapVirasoro.eval(b_p, τ) - BootstrapVirasoro.eval(b_m, τ)) /
            (2 * ϵ)
        block_der - byhand
    end

    byhand = (b_p(τ) - b_m(τ)) / (2 * ϵ)
    h = eval_series(b, τ)
    hprime_byhand = (eval_series(b_p, τ) - eval_series(b_m, τ)) / (2 * ϵ)

    @test isapprox(hprime_byhand, BootstrapVirasoro.eval_series_der(b, τ))

    @test b(τ) / eval_series(b, τ) ≈ q^(V0.dims[:left].δ) / etaDedekind(τ)

    @test BootstrapVirasoro.eval_der(b, τ) ≈
          q^(V0.dims.left.δ) / etaDedekind(τ) *
          (hprime_byhand + 2 * V0.dims.left.P * log(q) * h)

    test_vals =
        [big"0.3" + big"0.4" * im, big"0.6" + big"1.2" * im, big"10" + big"0.01" * im]
    for val in test_vals
        @test abs(eval_block_der(val)) < big"1e-23"
    end
end

@testset "Freg" begin
    V0 = Field(c, r = 2, s = 1)
    ϵ = 1e-40
    δ0ϵ = V0.dims.left.δ + ϵ
    V0ϵ = Field(c, diagonal = true, δ = δ0ϵ)
    V0_S2 = Field(c_S2, r = 2 * V0.r, s = V0.s)
    δ0_S2ϵ = V0_S2.dims.left.δ + ϵ
    V0_S2ϵ = Field(c_S2, diagonal = true, δ = δ0_S2ϵ)
    b = Block(co, :s, V0, Δmax)
    bϵ = Block(co, :s, V0ϵ, Δmax)
    b_S2 = Block(co_S2, :s, V0_S2, Δmax)
    b_S2ϵ = Block(co_S2, :s, V0ϵ, Δmax)

    lr_cache = BootstrapVirasoro.LRPosCache(τ, b)
    xlr_cache = BootstrapVirasoro.LRPosCache(x, b_S2)
    @test isapprox(
        bϵ.cblocks.left(τ),
        co.Rmn.left[V0.r, V0.s] / ϵ * BootstrapVirasoro.eval_lr_op(b, lr_cache)[1] +
        b.cblocks.left(τ),
        rtol = 1e-25,
    )# F_{Prs + ϵ} = R/ϵ F_{Pr,-s} + F^reg

    bϵ.cblocks.left(τ) - (
        co.Rmn.left[V0.r, V0.s] / ϵ * BootstrapVirasoro.eval_lr_op(b, lr_cache)[1] +
        b.cblocks.left(τ)
    )

    @test eval_series(b.cblocks.left, τ_cache.left) ≈
          eval_series(b_S2.cblocks.left, x)
    # regularized series

    Freg = b[:left](τ)
    Rrs = co.Rmn.left[V0.r, V0.s]
    Fminus = BootstrapVirasoro.eval_lr_op(b, lr_cache)[1]
    Freg_S2 = b_S2[:left](x)
    Rrs_S2 = co_S2.Rmn.left[:s][V0_S2.r, V0_S2.s]
    Fminus_S2 = BootstrapVirasoro.eval_lr_op(b_S2, xlr_cache)[1]

    @test Freg_S2 / prefac(b_S2[:left], x) ≈
          (Freg + 8log(big"2") * Rrs * Fminus) / prefac(b[:left], τ)
    # regularized block

    @test (Freg_S2 - 4log(big"2") * Rrs_S2 * Fminus_S2) / prefac(b_S2[:left], x) ≈
          Freg / prefac(b[:left], τ)
    # regularized block
end

@testset "Flog" begin
    V = Field(c, r = 2, s = 1)
    V_S2 = Field(c_S2, r = 2 * V.r, s = V.s)

    b = Block(co, :s, V, Δmax)
    b_op = Block(co, :s, Field(V.c, r = V.r, s = -V.s), Δmax)
    b_S2 = Block(co_S2, :s, V_S2, Δmax)
    b_op_S2 = Block(co_S2, :s, Field(V_S2.c, r = V_S2.r, s = -V_S2.s), Δmax)

    P23 = V.dims.left.P
    P43_S2 = V_S2.dims.left.P

    @test eval_series(b[:left], τ_cache.left) ≈ eval_series(b_S2[:left], x)
    # regularized series match

    missing_terms = b[:left].missing_terms
    missing_term = log(q) * Base.evalpoly(q, missing_terms)

    @test b[:left](τ) / total_prefactor(b[:left], τ) ≈
          eval_series(b[:left], τ_cache.left) + missing_term
    # the regularised block is as expected

    missing_terms_S2 = b_S2[:left].missing_terms
    missing_term_S2 = log(16q_S2) * Base.evalpoly(16q_S2, missing_terms_S2)

    @test 2 * Base.evalpoly(q, missing_terms) ≈
          Base.evalpoly(16q_S2, missing_terms_S2)

    @test log(q_S2) / log(16q_S2) * missing_term_S2 ≈ missing_term

    @test b_S2[:left](x) / total_prefactor(b_S2[:left], x) ≈
          eval_series(b_S2[:left], x) + missing_term_S2
    # the regularised block is as expected

    @test b[:left](τ) / total_prefactor(b[:left], τ) - missing_term ≈
          b_S2[:left](x) / total_prefactor(b_S2[:left], x) - missing_term_S2
end

setprecision(BigFloat, 40, base = 10)
c = CentralCharge(β = big"1.2" + big"0.1" * im)
c_S2 = CentralCharge(β = c.β / sqrt(big(2)))

chan_ind = [(2, 1), (4, 4), (3, 1), (5, 1), (5, 2)]
Vs = [Field(c, r = r, s = s) for (r, s) in chan_ind]
Vs_S2 = [Field(c_S2, r = 2 * V.r, s = V.s) for V in Vs]

ind = (0, big"0.3" + big"0.2" * im)
co = Correlation(Field(c, r = ind[1], s = ind[2]), Δmax)

ind_S2 = [(0, 1 // 2), (ind[1], ind[2] / 2), (0, 1 // 2), (0, 1 // 2)]
co_S2 = Correlation([Field(c_S2, r = r, s = s) for (r, s) in ind_S2], Δmax)

bs = [Block(co, V) for V in Vs]
bs_S2 = [Block(co_S2, :s, V_S2) for V_S2 in Vs_S2]

ell(b) = ell(b.corr.fields, b.chan_field.r, b.chan_field.s)

@testset "ell" begin
    # ell/sqrt(2) + 16 log(2) \beta'^{-1} s = ell^{S^2}
    for (i, b) in enumerate(bs)
        @test ell(b) / sqrt(big"2") + 16 * log(big"2") / c_S2.β * Vs[i].s ≈
              ell(bs_S2[i])
    end
end

prefs = [
    total_prefactor(b[:left], τ) * total_prefactor(b[:right], conj_q(τ, b.corr))
    for b in bs
]
prefs_S2 = [
    total_prefactor(b[:left], x) * total_prefactor(b[:right], conj_q(x, b.corr))
    for b in bs_S2
]

@testset "Full log block" begin
    for (i, b) in enumerate(bs)
        @test b(τ) / prefs[i] ≈ bs_S2[i](x) / prefs_S2[i]
    end
end

@testset "Interchiral" begin
    P = big"0.53" + big"0.11" * im
    V = Field(c, P = P, diagonal = true)
    V_S2 = Field(c_S2, P = sqrt(big"2") * P, diagonal = true)

    s = shift_D(co.fields, V)
    s_S2 = shift_D(co_S2.fields, V_S2)

    # shift(D^S2) = shift(D) * 16^(-8β^-1 (P - β^{-1}/2))
    @test isapprox(s_S2 / s / 16^(4 / c.β^2 * (V.s + 1)), 1)

    b = Block(co, :s, V, interchiral = true)
    prefactor =
        BootstrapVirasoro.PosCache(τ, b.blocks[1][:left]).prefactor *
        BootstrapVirasoro.PosCache(
            conj_q(τ, b.blocks[1].corr),
            b.blocks[1][:right],
        ).prefactor

    b_S2 = Block(co_S2, :s, V_S2, interchiral = true)
    prefactor_S2 =
        BootstrapVirasoro.PosCache(x, b_S2.blocks[1][:left]).prefactor *
        BootstrapVirasoro.PosCache(
            conj_q(x, b_S2.blocks[1].corr),
            b_S2.blocks[1][:right],
        ).prefactor

    for i = 1:3
        @test isapprox(
            b.blocks[i](τ) / prefactor,
            b_S2.blocks[i](x) / prefactor_S2 * 16^(-4b.fields[i].dims.left.P^2),
            rtol = 1e-20,
        )
    end

    @test isapprox(
        b(τ) / prefactor,
        b_S2(x) / prefactor_S2 * 16^(-4P^2),
        rtol = 1e-20,
    )
end

setprecision(BigFloat, 40, base = 10)
c = CentralCharge(β = big"1.2" + big"0.1" * im)
c_S2 = CentralCharge(β = c.β / sqrt(big(2)))
Δmax = 30

ind = (1, 0)
co = Correlation(Field(c, r = ind[1], s = ind[2]), Δmax)

ind_S2 = [(0, 1 // 2), (ind[1], ind[2] / 2), (0, 1 // 2), (0, 1 // 2)]
co_S2 = Correlation([Field(c_S2, r = r, s = s) for (r, s) in ind_S2], 2 * Δmax)

τ = big"0.36" + big"1.14" * im
x = xfromq(exp(im * (π * τ)))

xfromτ(τ) = xfromq(exp(im * (π * τ)))

function prefA(τ, co_S2, chan, lr)
    x_cache = PosCache(xfromτ(τ), co_S2[lr], chan)
    return inv(etaDedekind(τ) * x_cache.prefactor)
end

# non-chiral prefactor
ncprefA(τ, co_S2, chan) =
    prefA(τ, co_S2, chan, :left) * prefA(-conj(τ), co_S2, chan, :right)

@testset "Sphere <-> torus prefactors" begin
    for ind in [(5, 1), (11 // 2, 2 // 11), (3, 1), (2, 0), (1, 1), (3//2, 2//3)]
        V = Field(c, r = ind[1], s = ind[2])
        P, Pbar = V.dims.left.P, V.dims.right.P
        V_S2 = Field(c_S2, r = 2 * V.r, s = V.s)
        b = Block(co, V)
        b_S2 = @channels Block(co_S2, chan, V_S2)

        V1 = co.fields[1]
        Δ1, Δ1bar = V1.dims.left.Δ, V1.dims.right.Δ

        ds = co_S2[:left].fields
        e1 = -ds[1].Δ - ds[4].δ
        e2 = ds[1].Δ + sum(ds[i].δ for i = 2:4)
        dsbar = co_S2[:right].fields
        e1bar = -dsbar[1].Δ - dsbar[4].δ
        e2bar = dsbar[1].Δ + sum(dsbar[i].δ for i = 2:4)

        @test b(τ) * 16^(2P^2 + 2Pbar^2) ≈ b_S2.s(x) * ncprefA(τ, co_S2, :s)

        @test prefA(τ, co_S2, :s, :left) ≈
              prefA(-1 / τ, co_S2, :t, :left) * (-im * τ)^(-Δ1)

        @test prefA(τ, co_S2, :s, :left) ≈
              prefA((τ - 2) / (τ - 1), co_S2, :u, :left) *
              exp(im * (π * (e1 + e2))) *
              (-im * (τ - 1))^(-Δ1)

        @test prefA(-conj(τ), co_S2, :s, :right) ≈
              prefA(1 / conj(τ), co_S2, :t, :right) * (im * conj(τ))^(-Δ1bar)

        @test prefA(-conj(τ), co_S2, :s, :right) ≈
              prefA(-conj((τ - 2) / (τ - 1)), co_S2, :u, :right) *
              exp(-im * (π * (e1bar + e2bar))) *
              (im * (conj(τ) - 1))^(-Δ1bar)

        @test ncprefA(τ, co_S2, :s) ≈
              ncprefA(-1 / τ, co_S2, :t) * τ^-Δ1 * conj(τ)^-Δ1bar * im^spin(V1)

        @test ncprefA(τ, co_S2, :s) ≈
              ncprefA((τ - 2) / (τ - 1), co_S2, :u) *
              (τ - 1)^(-Δ1) *
              (conj(τ) - 1)^(-Δ1bar)

        f = ncprefA(τ, co_S2, :s)

        @test b(τ) * 16^(2P^2 + 2Pbar^2) ≈ b_S2.s(x) * f

        @test b(-1 / τ) * 16^(2P^2 + 2Pbar^2) ≈
              b_S2.t(1 - x) * f * τ^Δ1 * conj(τ)^Δ1bar * (-1)^spin(V1_S2)

        @test b(-1 / τ) * 16^(2P^2 + 2Pbar^2) ≈
              b_S2.t(1 - x) * ncprefA(-1 / τ, co_S2, :t)

        @test b((τ - 2) / (τ - 1)) * 16^(2P^2 + 2Pbar^2) ≈
              b_S2.u(1 / x) * f * (τ - 1)^Δ1 * (conj(τ) - 1)^Δ1bar

        @test b(τ) * 16^(2P^2 + 2Pbar^2) ≈ b_S2.s(x) * ncprefA(τ, co_S2, :s)

        @test b(-1 / τ) * 16^(2P^2 + 2Pbar^2) ≈
              b_S2.t(1 - x) * ncprefA(-1 / τ, co_S2, :t)

        @test b((τ - 2) / (τ - 1)) * 16^(2P^2 + 2Pbar^2) ≈
              b_S2.u(1 / x) * ncprefA((τ - 2) / (τ - 1), co_S2, :u)
    end
end

@testset "Interchiral" begin
    V = Field(c, r = 11 // 2, s = 2 // 11)
    P, Pbar = V.dims.left.P, V.dims.right.P
    V_S2 = Field(c_S2, r = 2 * V.r, s = V.s)
    b = Block(co, V, interchiral = true)
    b_S2 = @channels Block(co_S2, chan, V_S2, interchiral = true)
    chan = :u
    t = BootstrapVirasoro.modular_param(chan, τ)
    X = BootstrapVirasoro.crossratio(chan, x)
    @test b(t) * 16^(2P^2 + 2Pbar^2) ≈ b_S2[chan](X) * ncprefA(t, co_S2, chan)
    @test b.blocks[2](t) * 16^(2(P - 1 / c.β)^2 + 2(Pbar + 1 / c.β)^2) ≈
          b_S2[chan].blocks[2](X) * ncprefA(t, co_S2, chan)
    # shift(D^S2) = shift(D) * 16^(-8β^-1 (P - β^{-1}/2))
    @test b.shifts[2] * 16^(-4 / c.β^2 * (V.s + 1)) ≈ b_S2.u.shifts[2]

    for ind in [
        (5, 1),
        (5, 0),
        (4, 1 // 2),
        (9 // 2, 6 // 9),
        (11 // 2, 2 // 11),
        (3, 1),
        (2, 0),
        (1, 1),
    ]
        V = Field(c, r = ind[1], s = ind[2])
        P, Pbar = V.dims.left.P, V.dims.right.P
        V_S2 = Field(c_S2, r = 2 * V.r, s = V.s)
        b = Block(co, V, interchiral = true)
        b_S2 = @channels Block(co_S2, chan, V_S2, interchiral = true)
        for chan in BootstrapVirasoro.CHANNELS
            t = BootstrapVirasoro.modular_param(chan, τ)
            X = BootstrapVirasoro.crossratio(chan, x)
            @test isapprox(
                b(t) * 16^(2P^2 + 2Pbar^2),
                b_S2[chan](X) * ncprefA(t, co_S2, chan),
                rtol = 1e-23,
            )
        end
    end
end