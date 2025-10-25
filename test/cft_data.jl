@testset "Central Charges" begin
    #ensure the relation between b and β does not change
    c1 = CentralCharge(c = -1.1 + 0.2im)
    b = c1.b
    c2 = CentralCharge(b = b)
    @test c1.c == c2.c
    @test typeof(c1) == CC{ComplexF64}
    @test c1.β == c2.β

    c3 = CC(β = big"0.1")
    @test typeof(c3) == CC{Complex{BigFloat}}

    c3 = CC(β=1)
    @test typeof(c3) == CC{Complex{Float64}}
end

@testset "ConformalDimensions" begin
    c = CC(β = 1.0)
end

@testset "Fields" begin
    #ensure the relation between p and P does not change
    c1 = CentralCharge(c = -1.1 + 0.2im)
    V1 = Field(c1, P = 0.5)
    p = V1[:left].P
    V2 = Field(c1, P = p)
    @test all(isapprox.(V1[:left].P, V2[:left].P))

    #ensure the keyword diagonal also works for fields given from Kac indices
    V1 = Field(c1, r = 3, s = 4, diagonal = true)
    @test V1[:left].δ ≈ V1[:right].δ
    @test V1.degenerate

    V1 = Field(c1, r = 2, s = 5, diagonal = true)
    @test V1[:left].δ ≈ V1[:right].δ
    V2 = Field(c1, diagonal = true, r = 2, s = 5)
    @test V2 == V1
end
