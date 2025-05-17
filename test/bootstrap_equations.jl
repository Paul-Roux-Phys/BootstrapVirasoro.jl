using BootstrapVirasoro

c = CC(β = 1/(big"0.8" + big"0.1"*im))
V1 = Field(c, r=1//2, s=0)
Δmax = 15.
# setprecision(BigFloat, 13, base=10)
co = Correlation(V1, V1, V1, V1, Δmax=Δmax)
conditions(r, s) = r*s%1==0 && real(total_dimension(Field(c, r=r, s=s))) <= Δmax
indices = generate_pairs(1//2:1//2:30, conditions)
vs = [Field(c, r=r, s=s) for (r, s) in indices]
S = Dict(
    chan => Spectrum(co, chan, vs, interchiral=true, Δmax=Δmax)
    for chan in (:s, :t, :u)
)

xs = rand(Complex{BigFloat}, 10)

@time bootstrap_equations(S, xs)
@btime [bootstrap_equations(S, x) for x in xs]

V = Field(c, r=1, s=0)
b = Block(co, :s, V, interchiral=true, Δmax=Δmax)

x = big"0.3" + big"0.1"*im

c = CC(β = (big"0.8" + big"0.1"*im))
β = c.β
B = -c.β^(2)
function gamma_reg(x)
    """ Redefinition of the Gamma function, where the value
    at the pole is the residue. """
    if floor(real(x)) == x && real(x) <= 0
        return (-1)^floor(Int, real(x)) / gamma(1 - x)
    else
        return gamma(x)
    end
end

function MM(P, P1, P2)
    prod(
        gamma_reg(1 // 2 + β * (P + pm1 * P1 + pm2 * P2))
        for pm1 in (-1, 1)
        for pm2 in (-1, 1)
    )
end

function C_shift_test(V1, V2, V3)
    β = V1.c.β
    P1, P2, P3 = V1.P[:left], V2.P[:left], V3.P[:left]
    bP1, bP2, bP3 = V1.P[:right], V2.P[:right], V3.P[:right]
    prod(
        gamma_reg(1//2 - (bP3 + pm1 * bP2 + pm2 * bP1) / β) / 
            gamma_reg(1//2 - (P3 + pm1 * P2 + pm2 * P1) / β) 
        for pm1 in (-1, 1)
        for pm2 in (-1, 1)
    ) * β^(-4β^(-2)*V3.s)
end

import BootstrapVirasoro: shift_C123
V2 = Field(c, r=3//2, s=2//3)
shift_C123(V2, V2, V)
C_shift_test(V2, V2, V)

function D_ratio(r, s, V1, V2, V3, V4)
    P1, P2, P3, P4 = V1.P[:left], V2.P[:left], V3.P[:left], V4.P[:left]
    bP1, bP2, bP3, bP4 = V1.P[:right], V2.P[:right], V3.P[:right], V4.P[:right]
    prod(
        gamma_reg(-ϵ * B * r + η * s)^(-ϵ) *
        gamma_reg((1 - η) / 2 - (ϵ * r + η) * B - s)^(-ϵ)
        for ϵ in (-1, 1)
        for η in (-1, 1)
    ) * MM(P_rs(r, s, V1.c), P1, P2) / MM(-P_rs(r, s, V1.c), bP1, bP2) *
    MM(P_rs(r, s, V1.c), P3, P4) / MM(-P_rs(r, s, V1.c), bP3, bP4)
end
D_ratio(0, 1, co.fields...)
co.fields
c

S = (S_s, S_s, S_s)

size(S_s)
nb_blocks(S_s)

xs = BootstrapVirasoro.draw_points(10)

@btime bootstrap_equations(co, S, x)
BootstrapVirasoro.disable_distributed()
@time bootstrap_equations(co, S)
BootstrapVirasoro.enable_distributed()
@time bootstrap_equations(co, S)
M = bootstrap_equations(co, S, extrapoints=10)

M \ zeros(size(M)[1])

A = [[1, 2, 3], [4, 5, 6]]
Matrix(A)

length(S_s.blocks)

Vector{Vector{Float64}}(undef)

                      

for (i, j) in enumerate(S_s.blocks)
    println(i)
end


