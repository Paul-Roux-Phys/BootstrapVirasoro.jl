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
end

@testset "One point correlations" begin

end