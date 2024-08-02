const RmnTable{T} = Dict{Tuple{Int, Int}, Tuple{T, T}}
const CNmnTable{T} = Dict{Tuple{Int, Int, Int}, Tuple{T, T}}

double_prod_in_Dmn(m, n, B) = prod(prod((r^2*B - s^2/B)^2 for s in 1:n-1) for r in 1:m-1)

δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)

"""
Quantity `D_{m,n}` appearing in the denominator of the residues R_{m,n}
"""
function Dmn(m::Int, n::Int, B::Number)
    if m == 1 && n == 1 # treat cases m = 1, n=1 separately
        return 1
    elseif m == 1
        return n * prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
    elseif n == 1
        return m * prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
    else
        f1 = prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
        f2 = prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
        f3 = double_prod_in_Dmn(m, n, B)
        return m*n*f1*f2*f3
    end
end

include("four_point.jl")
include("one_point.jl")