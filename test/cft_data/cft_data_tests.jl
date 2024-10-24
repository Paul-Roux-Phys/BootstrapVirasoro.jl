@testset "Central Charges" begin
    #ensure the relation between b and β does not change
    c1 = CentralCharge(:c, -1.1 + 0.2im)
    b = c1.b
    c2 = CentralCharge(:b, b)
    @test c1.c == c2.c
    @test c1.β == c2.β
end

@testset "Fields" begin
    #ensure the relation between p and P does not change
    c1 = CentralCharge(:c, -1.1 + 0.2im)
    V1 = Field(c1, :P, 0.5, diagonal=true)
    p = V1.P[:left]
    V2 = Field(c1, :P, p, diagonal=true)
    @test V1.P == V2.P

    #ensure the keyword diagonal also works for fields given from Kac indices
    V1 = Field(c1, Kac=true, r=3, s=4, diagonal=true)
    @test V1.δ[:left] == V1.δ[:right]


    #ensure degenerate and diagonal work well together
    V1 = Field(c1, Kac=true, degenerate=true, r=2, s=5, diagonal=true)
    @test V1.δ[:left] == V1.δ[:right]
end