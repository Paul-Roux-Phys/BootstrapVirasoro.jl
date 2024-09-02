double_prod_in_Dmn(m, n, B) = prod(prod((r^2*B - s^2/B)^2 for s in 1:n-1) for r in 1:m-1)

function Dmn(m::Int, n::Int, B::Number)
    # treat cases m = 1, n=1 separately
    m == 1 && n == 1 && return 1
    m == 1 && return n * prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
    n == 1 && return m * prod(r^2 * B * (r^2 * B - n^2 / B) for r in 1:m-1)
    f1 = prod(r^2 * B * (r^2 * B - n^2 / B) for r in 1:m-1)
    f2 = prod(s^2 / B * (s^2 / B - m^2 * B) for s in 1:n-1)
    f3 = double_prod_in_Dmn(m, n, B)
    return m * n * f1 * f2 * f3
end

include("Residues/residues_1pt.jl")
include("Residues/residues_4pt.jl")

@memoize function computeCNmn(N, m, n, c, Rmn)
    B = c.B
    (!((m, n) in keys(Rmn)) || m*n > N) && return 0
    m*n == N && return Rmn[(m, n)]
    res = sum(sum(computeCNmn(N - m * n, mp, np, c, Rmn) / (δrs(m, -n, B) - δrs(mp, np, B))
                  for mp in 1:N-m*n if mp * np <= N - m * n)
              for np in 1:N-m*n)
    return Rmn[(m, n)] * res
end