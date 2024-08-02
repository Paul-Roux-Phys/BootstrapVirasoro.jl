@testset "Relation to sphere four-point blocks" begin
    c_torus = CentralCharge(:b, 1.2+.1*1im);
    c_sphere = CentralCharge(:b, (1.2+.1*1im)/sqrt(2))

    q = JuliVirBootstrap.FourPointBlocksSphere.qfromx(0.05)

    P = 0.23+.11im
    P1 = 0.41+1.03im
    V_torus_chan = Field(c_torus, :P, P, diagonal=true)
    δ_torus = V_torus_chan.δ[left]
    δ11_torus = Field(c_torus, Kac=true, r=1, s=1, diagonal=true).δ[left]
    V_torus_ext = Field(c_torus, :P, P1, diagonal=true)

    V_sphere_chan = Field(c_sphere, :P, sqrt(2)*P, diagonal=true)
    δ_sphere = V_sphere_chan.δ[left]
    δ21_sphere = Field(c_sphere, Kac=true, r=2, s=1, diagonal=true).δ[left]
    δ12_sphere = Field(c_sphere, Kac=true, r=1, s=2, diagonal=true).δ[left]
    V_sphere_ext = Field(c_sphere, :P, P1/sqrt(2), diagonal=true)
    VKac_sphere = Field(c_sphere, Kac=true, r=0, s=1//2, diagonal=true)

    corr_torus = OnePointCorrelation(c_torus, V_torus_ext)
    block_torus = OnePointBlockTorus(V_torus_chan)

    corr_sphere = FourPointCorrelation(c_sphere, [VKac_sphere, V_sphere_ext, VKac_sphere,VKac_sphere])
    block_sphere = FourPointBlockSphere(:s, V_sphere_chan)

    h1 = JuliVirBootstrap.OnePointBlocksTorus.H(q^2, 5, block_torus, corr_torus, left)
    h2 = JuliVirBootstrap.FourPointBlocksSphere.H(q, 5, block_sphere, corr_sphere, left)

    @test isapprox(h1, h2, atol=1e-12)

end