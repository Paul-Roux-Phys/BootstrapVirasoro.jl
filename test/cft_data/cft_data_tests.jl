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
    V1 = Field(c1, r=3, s=4, diagonal=true)
    @test V1.δ[:left] == V1.δ[:right]

    #ensure degenerate and diagonal work well together
    V1 = Field(c1, degenerate=true, r=2, s=5, diagonal=true)
    @test V1.δ[:left] == V1.δ[:right]
    # and independently
    V2 = Field(c1, degenerate=true, r=2, s=5)
    @test V2 == V1
    V3 = Field(c1, diagonal=true, r=2, s=5)
    @test V3 == V1

    # test shift()
    @test shift(V1, 2, :r).r == 4
    V1 = Field(c1, P=0.5, diagonal=true)
    @test shift(V1, 1).dim.P - V1.dim.P == 1/c1.β/2
end
