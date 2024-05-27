using JuliVirBootstrap
using Test

@testset "CFTData.jl" begin

    #ensure the relation between b and β does not change
    c1 = CentralCharge("c", -1.1+.2im)
    b = c1["b"]
    c2 = CentralCharge("b", b)
    @test c1["c"] == c2["c"]
    @test c1["β"] == c2["β"]

    #ensure the relation between p and P does not change
    left = 1
    right = 2
    V1 = Field(c1, "P", 0.5, diagonal=true)
    p = V1["p"][left]
    V2 = Field(c1, "p", p, diagonal=true)
    @test V1["P"] == V2["P"]

    #ensure the keyword diagonal also works for fields given from Kac indices
    V1 = Field(c1, Kac=true, r=3, s=4, diagonal=true)
    @test V1["Δ"][left] == V1["Δ"][right]


    #ensure degenerate and diagonal work well together
    V1 = Field(c1, Kac=true, degenerate=true, r=2, s=5, diagonal=true)
    @test V1["Δ"][left] == V1["Δ"][right]

end

@testset "FourPointCorrelationFunctions" begin

    left=1
    right=2

    c = CentralCharge("β", 1.2+.1*1im)
    V1 = Field(c, "Δ", 0.23+.11im, diagonal=true)
    V2 = Field(c, "Δ", 3.43, diagonal=true)
    V3 = Field(c, "Δ", 0.13, diagonal=true)
    V4 = Field(c, "Δ", 1.3, diagonal=true)
    corr = FourPointCorrelation(c, V1, V2, V3, V4)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.Rmn(2, 1, corr, "s", left),
                   0.31097697185245077-0.70523695127635733im, # value taken from Sylvain's code
                   atol=1e-8)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.Rmn(3, 3, corr, "t", left),
                   4.3964194233662846e-5-1.1534661157146291e-5im, # value taken from Sylvain's code
                   atol=1e-8)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.computeCNmn(7, 2, 3, corr, "s", left),
                   0.0019498393368877166+0.0026353877950837049im, # value taken from Sylvain's code
                   atol=1e-8)

end

@testset "FourPointBlocks" begin

    left=1;
    right=2;

    import JuliVirBootstrap.FourPointBlocksSphere.qfromx

    c_sphere = CentralCharge("b", (1.2+.1*1im)/sqrt(2))

    q = JuliVirBootstrap.FourPointBlocksSphere.qfromx(0.05)

    P = 0.23+.11im
    P1 = 0.41+1.03im

    V_sphere_chan = Field(c_sphere, "P", sqrt(2)*P, diagonal=true)
    V_sphere_ext = Field(c_sphere, "P", P1/sqrt(2), diagonal=true)
    VKac_sphere = Field(c_sphere, Kac=true, r=0, s=1//2, diagonal=true)

    corr_sphere = FourPointCorrelation(c_sphere, [VKac_sphere, V_sphere_ext, VKac_sphere,VKac_sphere])
    block_sphere = FourPointBlockSphere("s", V_sphere_chan)

    h = JuliVirBootstrap.FourPointBlocksSphere.H(q, 5, block_sphere, corr_sphere, left)

    @test isapprox(h, 0.9999955375834808 - 2.735498726466085e-6im, atol=1e-8) # value from Sylvain's code


    setprecision(BigFloat, 64)

    c = CentralCharge("β", big(1.2+.1*1im));
    V1 = Field(c, "Δ", 0.23+.11im, diagonal=true);
    V2 = Field(c, "Δ", 3.43, diagonal=true);
    V3 = Field(c, "Δ", 0.13, diagonal=true);
    V4 = Field(c, "Δ", 1.3, diagonal=true);
    V = Field(c, "Δ", 0.1, diagonal = true);

    corr = FourPointCorrelation(c, [V1, V2, V3, V4])

    bl_s = FourPointBlockSphere("s", V)
    bl_t = FourPointBlockSphere("t", V)
    bl_u = FourPointBlockSphere("u", V)

    x=0.05

    # comparing to values from Sylvain's code
    @test isapprox(JuliVirBootstrap.FourPointBlocksSphere.block_chiral(x, 6, bl_s, corr, left), 2337.4038141240320199350204984981259378760811288542 + 4771.3912725970751669197262259253749217475400016186im, rtol = 1e-10)
    @test isapprox(JuliVirBootstrap.FourPointBlocksSphere.block_chiral(x, 6, bl_t, corr, left), 52191.790807047848992452669811987274395806031692488 - 140430.98553278617162374003412214159828722759436549im,rtol = 1e-10)
    @test isapprox(JuliVirBootstrap.FourPointBlocksSphere.block_chiral(x, 6, bl_u, corr, left), 852.92814340196565010929995606986011067184449511918 + 359.96303529282323934093142050535102602840290239155im, rtol = 1e-10)

end

@testset "OnePointBlocks" begin
    left=1;
    right=2;

    import JuliVirBootstrap.FourPointBlocksSphere.qfromx
    c_torus = CentralCharge("b", 1.2+.1*1im);
    c_sphere = CentralCharge("b", (1.2+.1*1im)/sqrt(2))

    q = JuliVirBootstrap.FourPointBlocksSphere.qfromx(0.05)

    P = 0.23+.11im
    P1 = 0.41+1.03im
    V_torus_chan = Field(c_torus, "P", P, diagonal=true)
    δ_torus = V_torus_chan["δ"][left]
    δ11_torus = Field(c_torus, Kac=true, r=1, s=1, diagonal=true)["δ"][left]
    V_torus_ext = Field(c_torus, "P", P1, diagonal=true)

    V_sphere_chan = Field(c_sphere, "P", sqrt(2)*P, diagonal=true)
    δ_sphere = V_sphere_chan["δ"][left]
    δ21_sphere = Field(c_sphere, Kac=true, r=2, s=1, diagonal=true)["δ"][left]
    δ12_sphere = Field(c_sphere, Kac=true, r=1, s=2, diagonal=true)["δ"][left]
    V_sphere_ext = Field(c_sphere, "P", P1/sqrt(2), diagonal=true)
    VKac_sphere = Field(c_sphere, Kac=true, r=0, s=1//2, diagonal=true)

    corr_torus = OnePointCorrelation(c_torus, V_torus_ext)
    block_torus = OnePointBlockTorus(V_torus_chan)

    corr_sphere = FourPointCorrelation(c_sphere, [VKac_sphere, V_sphere_ext, VKac_sphere,VKac_sphere])
    block_sphere = FourPointBlockSphere("s", V_sphere_chan)

    h1 = JuliVirBootstrap.OnePointBlocksTorus.H(q^2, 5, block_torus, corr_torus, left)
    h2 = JuliVirBootstrap.FourPointBlocksSphere.H(q, 5, block_sphere, corr_sphere, left)

    @test isapprox(h1, h2, atol=1e-12)
end
