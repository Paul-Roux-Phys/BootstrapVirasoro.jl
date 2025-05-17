@testset "Four point correlations" begin
    @testset "Normal residues" begin
        c = CentralCharge(:c, 0.5)
        V1 = Field(c, r=1, s=0)
        V2 = Field(c, r=2, s=0)
        corr = Correlation(V1, V2, V2, V1, 6)

        # litteral values are taken from Sylvain's code
        @test isapprox(corr._Rmn[:left][:s][(3, 2)], -3.11111e-7,
            atol=1e-8)

        c = CentralCharge(:c, big"0.1")
        V1 = Field(c, r=1 // 2, s=0)
        V2 = Field(c, r=3 // 2, s=0)
        corr = Correlation(V1, V2, V2, V1, 10)

        @test isapprox(
            corr._CNmn[:left][:t][(4, 1, 4)],
            big"-0.00004010759186977254114462018639739857292228437487",
            atol=1e-20
        )
    end

    @testset "Regularised residues as limits of normal residues" begin
        import BootstrapVirasoro: computeRmn,
            Rmn_term,
            Rmn_zero_order,
            Rmn_term_vanishes,
            Dmn
        setprecision(BigFloat, 256)
        c = CentralCharge(:β, big"0.8" + big"0.1"*im)
        ϵ = 1 // big"10"^20
        
        V1 = Field(c, r=0, s=1);
        V2 = Field(c, r=0, s=1//2);
        V3 = Field(c, r=2, s=1//2);
        V4(ϵ) = Field(c, r=2, s=3//2 + ϵ);
        Vs(ϵ) = (V1, V2, V3, V4(ϵ))
        
        dls(ϵ) = Tuple(v.dims[:left] for v in Vs(ϵ))
        drs(ϵ) = Tuple(v.dims[:right] for v in Vs(ϵ))
        
        redirect_stderr(devnull) do
            Rreg_1_2 = computeRmn(1, 2, dls(0))

            P3, P4 = dls(ϵ)[3].P, dls(ϵ)[4].P
            P(r, s) = ConformalDimension(c, r=r, s=s).P
            Rϵ12 = computeRmn(1, 2, dls(ϵ)) / (P3 - P4 + P(0, 1))
            @test isapprox(Rreg_1_2, Rϵ12, rtol=1e-18)
        end

        R12(ϵ) = computeRmn(1, 12, dls(ϵ))
        barR12(ϵ) = computeRmn(1, 12, drs(ϵ))

        redirect_stderr(devnull) do 
            Rratio_reg = (-1)^Rmn_zero_order(1, 12, dls(0)) * R12(0) / barR12(0)
            Rratio_ϵ = R12(ϵ) / barR12(ϵ)

            @test isapprox(Rratio_reg, Rratio_ϵ, rtol=1e-18)
        end
    end
end
