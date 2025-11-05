using BootstrapVirasoro.LoopModels
import BootstrapVirasoro: xfromq, conj_q, PosCache, etaDedekind
import BootstrapVirasoro.LoopModels: shift_D
setprecision(BigFloat, 40, base=10)

@testset "ConformalDimensions" begin
    c = CC(β = 0.8+0.3im)
    d = CD(c, P = 0.8+0.3im)
    @test shift(d, 2).P - d.P ≈ -1/c.β # s = -2βP
end

@testset "ConformalDimensions" begin
    c = CC(β = 0.8+0.3im)
    d = CD(c, P = 0.8+0.3im)
    @test shift(d, 2).s - d.s == 2
    @test shift(d, 2).P - d.P ≈ -1/c.β # s = -2βP
end

@testset "Fields" begin
    #ensure the relation between p and P does not change
    c = CC(c = -1.1 + 0.2im)
    V = Field(c, r = 2, s = 5, diagonal = true)
    # test shift()
    @test shift(V, 2).r == 2
    @test shift(V, 2)[:right].s == 7
    V = Field(c, P = 0.5, diagonal = true)
    @test abs(shift(V, 1)[:left].P - V[:left].P + 1/c.β/2) < 1e-14
    V = Field(c, r=2, s=5)
    @test shift(V, 2)[:right].s == -7
end

@testset "Interchiral 4pt" begin
        c = CC(β=-big"0.8" - big"0.1" * im)
        field1 = Field(c, r=1 // 2, s=0)
        field2 = Field(c, r=1, s=0)
        fields = [field2, field2, field1, field1]
        Δmax = 40
        co = Correlation(fields..., Δmax)
        x = big"0.4" + big"0.2" * im

        # interchiral, logarithmic
        J = Field(c, r=1, s=1)
        b = IBlock(co, :s, J, Δmax)
        @test b(x) ≈
              big"1.320255511354332911164464364785819105156" +
              big"0.4186925417664498703779197115719248226952" * im

        # interchiral, non-diagonal
        V3 = Field(c, r=3, s=1 // 3)
        b = IBlock(co, :s, V3, Δmax)
        @test b(x) ≈
              big"0.2052943176316875457496459173129386291016" -
              big"0.2003078699572151816572767384428219647201" * im

        # interchiral, degenerate
        id = Field(c, r=1, s=1, diagonal=true)
        b = IBlock(co, :s, id, Δmax)
        @test b(x) ≈
              big"1.439312717815166500340922134926376204051" +
              big"0.4561207475025284508099330175652862628828" * im

        # interchiral, generic diagonal
        V4 = Field(c, r=0, s=big"0.5" + big"0.3" * im)
        b = IBlock(co, :s, V4, Δmax)
        @test b(x) ≈
              big"1.396132665254477154107379890952326016033" +
              big"0.184273048386930095042991719005258073884" * im

        V5 = Field(c, r=2, s=0)
        b = IBlock(co, :s, V5, Δmax)
        @test b(x) ≈
              big"0.7505040963332944454258635005057715597006" +
              big"0.04344655018615009393069583429415501260266" * im
end

const Sig = Channels{Rational}

c = CC(β=1 / (big"0.8" + big"0.1" * im))
ndiag_indices = [
        (r, s) for r = (1//2):(1//2):20 for s = (-1+1//(2r)):(1//(2r)):1 if r * s % 1 == 0
]
diag_field = Field(c, r=0, s=big"0.4" + big"0.1" * im)
fields = vcat([Field(c, r=r, s=s) for (r, s) in ndiag_indices], diag_field)

# Filter keeps only the fields/blocks that fulfill the condition given in the first argument.
function chan_parities(co::Correlation4)
        "Determine the parity of the number of legs in 4pt channels"
        V1, V2, V3, V4 = co.fields
        chan_parities =
                Channels{Rational}((V1.r + V2.r) % 1, (V1.r + V4.r) % 1, (V1.r + V3.r) % 1)
end

function LoopSpectra(co, fields, fs; parity=0)
        Vs = @channels filter(V -> V.r % 1 == chan_parities(co)[chan], fields)
        @channels ChannelSpectrum(co, chan, Vs[chan], fs[chan])
end

function precompute_blocks(co, fields; parity)
        parity != 0 &&
                (fields = filter(V -> V.diagonal || 0 <= V.s * abs(parity), fields))
        LoopSpectra(
                co,
                fields,
                Channels(
                        chan -> (V -> IBlock(co, chan, V, parity=parity)),
                ),
        )
end

function solve(specs, signature)
        specs = @channels filter(V -> V.r >= signature[chan], specs[chan])
        sys = BootstrapSystem(specs)
        evaluate_blocks!(sys)
        compute_linear_system!(sys)
        solve!(sys)
        return sys
end

ind = ((1 // 2, 0), (1 // 2, 0), (1 // 2, 0), (1 // 2, 0))
Δmax = 20
setprecision(BigFloat, 15, base=10)
co = Correlation([Field(c, r=r, s=s) for (r, s) in ind], Δmax)
blocks_even = precompute_blocks(co, fields, parity=1)

sig = Channels{Rational}(0, 1, 1)
sol = solve(blocks_even, sig)

@testset "Crossing" begin
        for chan in (:s, :t, :u)
                for (V, err) in sol.str_cst[chan].errors
                        @test abs(err * sol.str_cst[chan].constants[V]) < 1e-4
                end
        end
end

setprecision(BigFloat, 40, base=10)
c = CentralCharge(β=big"1.2" + big"0.1" * im)
c_S2 = CentralCharge(β=c.β / sqrt(big(2)))
Δmax = 30

ind = (1, 0)
co = Correlation(Field(c, r=ind[1], s=ind[2]), Δmax)

ind_S2 = [(0, 1 // 2), (ind[1], ind[2] / 2), (0, 1 // 2), (0, 1 // 2)]
co_S2 = Correlation([Field(c_S2, r=r, s=s) for (r, s) in ind_S2], 2 * Δmax)

τ = big"0.36" + big"1.14" * im
x = xfromq(exp(im * (π * τ)))

@testset "Interchiral 1pt" begin
        P = big"0.53" + big"0.11" * im
        V = Field(c, P=P, diagonal=true)
        V_S2 = Field(c_S2, P=sqrt(big"2") * P, diagonal=true)

        s = shift_D(co.fields, V)
        s_S2 = shift_D(co_S2.fields, V_S2)

        # shift(D^S2) = shift(D) * 16^(-8β^-1 (P - β^{-1}/2))
        @test isapprox(s_S2 / s / 16^(4 / c.β^2 * (V.s + 1)), 1)

        b = IBlock(co, :s, V)
        prefactor =
                BootstrapVirasoro.PosCache(τ, b.blocks[1][:left]).prefactor *
                BootstrapVirasoro.PosCache(conj_q(τ, b.blocks[1].corr), b.blocks[1][:right]).prefactor

        b_S2 = IBlock(co_S2, :s, V_S2)
        prefactor_S2 =
                BootstrapVirasoro.PosCache(x, b_S2.blocks[1][:left]).prefactor *
                BootstrapVirasoro.PosCache(
                        conj_q(x, b_S2.blocks[1].corr),
                        b_S2.blocks[1][:right],
                ).prefactor

        for i = 1:3
                @test isapprox(
                        b.blocks[i](τ) / prefactor,
                        b_S2.blocks[i](x) / prefactor_S2 *
                        16^(-4(b.blocks[i].chan_field[:left].P)^2),
                        rtol=1e-20,
                )
        end

        @test isapprox(
                b(τ) / prefactor,
                b_S2(x) / prefactor_S2 * 16^(-4P^2),
                rtol=1e-20,
        )
end

function prefA(τ, co_S2, chan, lr)
        x_cache = PosCache(xfromτ(τ), co_S2[lr], chan)
        return inv(etaDedekind(τ) * x_cache.prefactor)
end

# non-chiral prefactor
ncprefA(τ, co_S2, chan) =
        prefA(τ, co_S2, chan, :left) * prefA(-conj(τ), co_S2, chan, :right)

@testset "Interchiral" begin
        V = Field(c, r=11 // 2, s=2 // 11)
        P, Pbar = V[:left].P, V[:right].P
        V_S2 = Field(c_S2, r=2 * V.r, s=V.s)
        b = IBlock(co, V)
        b_S2 = @channels IBlock(co_S2, chan, V_S2)
        chan = :u
        t = BootstrapVirasoro.modular_param(chan, τ)
        X = BootstrapVirasoro.crossratio(chan, x)
        @test b(t) * 16^(2P^2 + 2Pbar^2) ≈ b_S2[chan](X) * ncprefA(t, co_S2, chan)
        @test b.blocks[2](t) * 16^(2(P - 1 / c.β)^2 + 2(Pbar + 1 / c.β)^2) ≈
              b_S2[chan].blocks[2](X) * ncprefA(t, co_S2, chan)
        # shift(D^S2) = shift(D) * 16^(-8β^-1 (P - β^{-1}/2))
        @test b.coeffs[2] * 16^(-4 / c.β^2 * (V.s + 1)) ≈ b_S2.u.coeffs[2]

        for ind in [
                (5, 1),
                (5, 0),
                (4, 1 // 2),
                (9 // 2, 6 // 9),
                (11 // 2, 2 // 11),
                (3, 1),
                (2, 0),
                (1, 1),
                (1, 1//2),
                (1, -1//2),
                (1//2, 1),
                (3//2, 1)
        ]
                V = Field(c, r=ind[1], s=ind[2])
                P, Pbar = V[:left].P, V[:right].P
                V_S2 = Field(c_S2, r=2 * V.r, s=V.s)
                b = IBlock(co, V)
                b_S2 = @channels IBlock(co_S2, chan, V_S2)
                for chan in BootstrapVirasoro.CHANNELS
                        t = BootstrapVirasoro.modular_param(chan, τ)
                        X = BootstrapVirasoro.crossratio(chan, x)
                        @test isapprox(
                                b(t) * 16^(2P^2 + 2Pbar^2),
                                b_S2[chan](X) * ncprefA(t, co_S2, chan) *
                                (chan === :u ? (-1)^(spin(V_S2)) : 1), 
                                rtol=1e-23,
                        )
                end
        end
end
