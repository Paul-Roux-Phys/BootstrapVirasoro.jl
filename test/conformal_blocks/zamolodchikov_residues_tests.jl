@testset "Four point correlations" begin
    c = CentralCharge(:c, 0.5)
    V1 = Field(c, Kac=true, r=1, s=0)
    V2 = Field(c, Kac=true, r=2, s=0)
    corr = Correlation(V1, V2, V2, V1, 6)

    # litteral values are taken from Sylvain's code
    @test isapprox(corr._Rmn[:left][:s][(3, 2)], -3.11111e-7,
        atol=1e-8)

    c = CentralCharge(:c, big"0.1")
    V1 = Field(c, Kac=true, r=1 // 2, s=0)
    V2 = Field(c, Kac=true, r=3 // 2, s=0)
    corr = Correlation(V1, V2, V2, V1, 10)

    @test isapprox(
        corr._CNmn[:left][:t][(4, 1, 4)],
        big"-0.00004010759186977254114462018639739857292228437487",
        atol=1e-20
    )

    # @testset "Regularised residues as limits of normal residues" begin
    #     c = CentralCharge(:β, big"0.8"+big"0.1"*im)
    #     Nmax = 15
    #     V1 = Field(c, diagonal=true, :Δ, 0.43)
    #     V2 = Field(c, Kac=true, r=3, s=2 // 3)
    #     V3 = Field(c, Kac=true, r=2, s=1 // 2)
    #     V4 = Field(c, Kac=true, r=2, s=3 // 2)

    #     V = Field(c, Kac=true, r=3, s=3)

    #     co = Correlation(V1, V2, V3, V4, Nmax)
    #     R1 = co._Rmn_reg[:left][:s][(1, 2)]
    #     R2 = co._Rmn_reg[:right][:s][(1, 2)]

        function shift_indices(V::FourFields, ϵ)
            s_shift = zeros(Rational, 4)
            for i in 1:4
                if i % 2 == 0
                    s_shift[i] = ϵ
                end
            end

            return [
                Field(co.c, Kac=true, r=V[i].r, s=V[i].s + s_shift[i])
                for i in 1:4
            ]
        end

        @test begin
            res = true
            for pair in keys(co._Rmn[:left][:s])
                if !(haskey(co._Rmn[:left][:s], pair))
                    res = false
                end
            end
            res
        end

        ϵ = big"1" // big"10"^10
        co_shift = Correlation(shift_indices((V1, V2, V3, V4), ϵ)..., Nmax)
        R3 = co_shift._Rmn[:left][:s][(1, 2)]
        R4 = co_shift._Rmn[:right][:s][(1, 2)]
    # end
end

@testset "One point correlations" begin

end