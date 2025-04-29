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

function evalpoly_buf!(buf_mul, buf, q, a)
    add_to!!(buf, Complex{BigFloat}(0), a[end]) # result = a_n
    for i in length(a)-1:-1:1 # result = result * q + a[n-1]
        mul_to_buf!!(buf_mul, buf, q)
        add_to!!(buf, buf, a[i])
    end
    deepcopy(buf)
end

function evalpoly_buf(q::Complex{BigFloat}, a)
    length(a) == 0 && return 0
    buf_mul = typeof(q)(0)
    buf = typeof(q)(0)
    evalpoly_buf!(buf_mul, buf, q, a)
end

function evalpoly_buf(q::Union{ComplexF64, Float64}, a)
    evalpoly(q, a)
end

function evapoly_buf(q::BigFloat, a)
    evalpoly(complex(q), q)
end

function evalpoly_buf(qs::Vector{T}, a::Vector{T}) where {T}
    nt = Threads.nthreads()
    buf_muls = [typeof(qs[1])(0) for _ in 1:nt]
    bufs = [typeof(qs[1])(0) for _ in 1:nt]
    res = Vector{typeof(qs[1])}(undef, length(qs))
    Threads.@threads for i in eachindex(qs)
        tid = Threads.threadid()
        res[i] = evalpoly_buf!(buf_muls[tid], bufs[tid], qs[i], a)
    end
    res
end

function evalpoly_buf(qs::Vector{T}, as::Vector{Vector{T}}) where {T}
    if use_distributed()
        pmap(a -> evalpoly_buf(qs, a), as)
    else
        [evalpoly_buf(qs, a) for a in as]
    end
end
