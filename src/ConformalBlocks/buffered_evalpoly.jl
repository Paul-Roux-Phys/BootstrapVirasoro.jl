import MutableArithmetics: mul_to!!, add_to!!, sub!!, broadcast!!, copy_if_mutable, add!!, sub_to!!

const BUF_MULS = [Complex{BigFloat}(0) for _ in Threads.nthreads()]
const BUFS = [Complex{BigFloat}(0) for _ in Threads.nthreads()]

#=
Computing series normally allocates a new number for each new term.
When using arbitrary precision numbers, this quickly grows to Gb of
memory.
We use in-place arithmetics to prevent this.
Complex multiplication requires an intermediary buffer, when computing
many products we use the same buffer throughout when possible.
Typically we divide memory usage by 3-4 orders of magnitude.
=#

function add_to!!(a::Complex, b::Complex, c::Complex)
    add_to!!(a.re, b.re, c.re)
    add_to!!(a.im, b.im, c.im)
end

# a *= q
function mul_to_buf!!(buf::Complex, a::Complex, q::Complex)
    mul_to!!(buf.re, a.re, q.re)
    mul_to!!(buf.im, a.re, q.im)
    mul_to!!(a.re, a.im, q.im)
    sub_to!!(a.re, buf.re, a.re)
    mul_to!!(a.im, a.im, q.re)
    add!!(a.im, buf.im)
    a
end

# --------------------
# Core in-place evaluation using Horner's method
# --------------------
function evalpoly_buf!(buf_mul, buf, q, a)
    add_to!!(buf, zero(eltype(a)), a[end])  # buf = a_n
    for i in length(a)-1:-1:1
        mul_to_buf!!(buf_mul, buf, q)       # buf *= q
        add_to!!(buf, buf, a[i])            # buf += a[i]
    end
    return deepcopy(buf)
end

# --------------------
# Scalar evaluation
# --------------------
function evalpoly_buf(q, a::AbstractVector)
    isempty(a) && return zero(promote_type(typeof(q), eltype(a)))
    T = promote_type(typeof(q), eltype(a))
    qT = convert(T, q)
    aT = convert(Vector{T}, a)
    buf_mul = zero(T)
    buf = zero(T)
    return evalpoly_buf!(buf_mul, buf, qT, aT)
end

# --------------------
# Fast path for built-in real/complex types (Base.evalpoly)
# --------------------
function evalpoly_buf(q::Union{Float64, ComplexF64}, a::AbstractVector{<:Union{Float64, ComplexF64}})
    evalpoly(q, a)
end

# --------------------
# Batched evaluation at multiple points
# --------------------
function evalpoly_buf(qs::Vector, a::Vector)
    isempty(qs) && return Vector{promote_type(eltype(qs), eltype(a))}()
    T = promote_type(eltype(qs), eltype(a))
    qsT = convert(Vector{T}, qs)
    aT = convert(Vector{T}, a)
    nt = Threads.nthreads()
    buf_muls = [zero(T) for _ in 1:nt]
    bufs = [zero(T) for _ in 1:nt]
    res = Vector{T}(undef, length(qs))
    Threads.@threads for i in eachindex(qsT)
        tid = Threads.threadid()
        res[i] = evalpoly_buf!(buf_muls[tid], bufs[tid], qsT[i], aT)
    end
    return res
end

# --------------------
# Batched evaluation of multiple polynomials
# --------------------
function evalpoly_buf(qs::Vector, polys::Vector{<:AbstractVector})
    isempty(polys) && return []
    T = promote_type(eltype(qs), eltype(first(polys)))
    qsT = convert(Vector{T}, qs)
    polysT = [convert(Vector{T}, a) for a in polys]
    # if use_distributed()
    #     return pmap(a -> evalpoly_buf(qsT, a), polysT)
    # else
        return [evalpoly_buf(qsT, a) for a in polysT]
    # end
end
