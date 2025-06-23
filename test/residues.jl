@testset "Normal residues" begin
    c = CentralCharge(c = 0.5)
    V1 = Field(c, r = 1, s = 0)
    V2 = Field(c, r = 2, s = 0)
    corr = Correlation(V1, V2, V2, V1, 6)

    # litteral values are taken from Sylvain's code
    @test isapprox(corr._Rmn[:left][:s][3, 2], -3.11111e-7, atol = 1e-8)

    c = CentralCharge(:c, big"0.1")
    V1 = Field(c, r = 1 // 2, s = 0)
    V2 = Field(c, r = 3 // 2, s = 0)
    corr = Correlation(V1, V2, V2, V1, 10)

    @test isapprox(
        corr._CNmn[:left][:t][4, 1, 4],
        big"-0.00004010759186977254114462018639739857292228437487",
        atol = 1e-20,
    )
end

@testset "Regularised residues as limits of residues" begin
    import BootstrapVirasoro: computeRmns, Rmn_term, Rmn_zero_order, Rmn_term_vanishes
    setprecision(BigFloat, 256)
    c = CentralCharge(:β, big"0.8" + big"0.1" * im)
    ϵ = 1 // big"10"^20

    V1 = Field(c, r = 0, s = 1)
    V2 = Field(c, r = 0, s = 1 // 2)
    V3 = Field(c, r = 2, s = 1 // 2)
    V4(ϵ) = Field(c, r = 2, s = 3 // 2 + ϵ)
    Vs(ϵ) = (V1, V2, V3, V4(ϵ))

    dls(ϵ) = Tuple(v.dims[:left] for v in Vs(ϵ))
    drs(ϵ) = Tuple(v.dims[:right] for v in Vs(ϵ))

    redirect_stderr(devnull) do
        P3, P4 = dls(ϵ)[3].P, dls(ϵ)[4].P
        P(r, s) = P_rs(r, s, c)

        Nmax = 12
        T = Complex{BigFloat}
        DRs = Matrix{T}(undef, (Nmax, Nmax))
        Pns = Matrix{T}(undef, (Nmax, Nmax))
        factors = Matrix{T}(undef, (Nmax, 2Nmax))

        Rl_reg = computeRmns(Nmax, dls(0))[2]
        Rlϵ = computeRmns(Nmax, dls(ϵ))[1]
        Rr_reg = computeRmns(Nmax, drs(0))[2]
        Rrϵ = computeRmns(Nmax, drs(ϵ))[1]

        Rlreg_1_2 = Rl_reg[1, 2]
        Rlϵ_1_2 = Rlϵ[1, 2] / (P3 - P4 + P(0, 1))

        @test isapprox(Rlreg_1_2, Rlϵ_1_2, rtol = 1e-18)

        Rreg_12 = Rl_reg[1, 12]
        barRreg_12 = Rr_reg[1, 12]
        Rϵ_12 = Rlϵ[1, 12]
        barRϵ_12 = Rrϵ[1, 12] # / (P3 - P4 + P(0, 1))

        Rratio_reg = (-1)^Rmn_zero_order(1, 12, dls(0)) * Rreg_12 / barRreg_12
        Rratio_ϵ = Rϵ_12 / barRϵ_12

        @test isapprox(Rratio_reg, Rratio_ϵ, rtol = 1e-18)
    end
end
