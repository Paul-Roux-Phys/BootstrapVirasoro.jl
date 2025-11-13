#=
# Allocating arbitrary precision numbers is expensive.
# Since julia number types are all immutable, whenever we do an operation like
# c = a + b
# Julia allocates a new number c even if c was previously allocated.
# For the performance critical parts of the code, i.e. computing blocks coefficients
# and evaluating block series, we use in-place arithmetics via the MutableArithmetics package API,
# if the types are arbitrary precision.
# We use global thread-local buffers to store results of intermediary computations, and implement
# the MutableArithmetics API for the Complex{BigFloat} numbers.
# Since the resulting syntax is not very readable, we only use this in performance critical parts
# of the library.
=#


using MutableArithmetics
const MA = MutableArithmetics
import MutableArithmetics: mutability, mutable_copy, buffer_for,
    promote_operation, operate_to!, operate!

const cplx_buffer  = Ref{Vector{Complex{BigFloat}}}()
const cplx_buffer2 = Ref{Vector{Complex{BigFloat}}}()
const real_buffer  = Ref{Vector{BigFloat}}()
const real_buffer2 = Ref{Vector{BigFloat}}()

# before Julia 1.8, there was only one global thread pool, in newer julia versions
# the threadid can be > nthreads.
@static if VERSION >= v"1.9"
    max_thread_id() = Threads.maxthreadid()
else
    max_thread_id() = Threads.nthreads()
end

function __init__()
    cplx_buffer[]  = [zero(Complex{BigFloat}) for _ in 1:max_thread_id()]
    cplx_buffer2[] = [zero(Complex{BigFloat}) for _ in 1:max_thread_id()]
    real_buffer[]  = [zero(BigFloat) for _ in 1:max_thread_id()]
    real_buffer2[] = [zero(BigFloat) for _ in 1:max_thread_id()]
end

MA.mutability(::Type{Complex{BigFloat}}) = IsMutable()
MA.mutable_copy(x::Complex{BigFloat}) = x

# copy
MA.promote_operation(::typeof(copy), ::Type{Complex{BigFloat}}) = Complex{BigFloat}

function MA.operate_to!(out::Complex{BigFloat}, ::typeof(copy), in::Complex{BigFloat})
    operate_to!(out.re, copy, in.re)
    operate_to!(out.im, copy, in.im)
    return out
end

MA.operate!(::typeof(copy), z::Complex{BigFloat}) = z 

# zero / one

MA.promote_operation(::typeof(zero), ::Type{Complex{BigFloat}}) = Complex{BigFloat}

MA.promote_operation(::typeof(one), ::Type{Complex{BigFloat}}) = Complex{BigFloat}

function MA.operate!(::typeof(zero), z::Complex{BigFloat})
    operate!(zero, z.re)
    operate!(zero, z.im)
    return z
end

function MA.operate!(::typeof(one), z::Complex{BigFloat})
    operate!(zero, z.re)
    operate!(zero, z.im)  # complex one should be (1,0)
    operate!(one, z.re)
    return z
end

# addition

MA.promote_operation(::typeof(+), ::Type{Complex{BigFloat}}, ::Type{Complex{BigFloat}}) = Complex{BigFloat}

function MA.operate_to!(out::Complex{BigFloat}, ::typeof(+), a::Complex{BigFloat}, b::Complex{BigFloat})
    operate_to!(out.re, +, a.re, b.re)
    operate_to!(out.im, +, a.im, b.im)
    return out
end

MA.operate_to!(out::Complex{BigFloat}, ::typeof(+), a::Complex{BigFloat}) = operate_to!(out, copy, a)
MA.operate!(::typeof(+), x::Complex{BigFloat}, y::Complex{BigFloat}) = operate_to!(x, +, x, y)
MA.operate!(::typeof(+), x::Complex{BigFloat}) = x  # unary + no-op

# subtraction
MA.promote_operation(::typeof(-), ::Vararg{Type{Complex{BigFloat}},N}) where {N} = Complex{BigFloat}

function MA.operate_to!(out::Complex{BigFloat}, ::typeof(-), a::Complex{BigFloat}, b::Complex{BigFloat})
    operate_to!(out.re, -, a.re, b.re)
    operate_to!(out.im, -, a.im, b.im)
    return out
end

MA.operate!(::typeof(-), z::Complex{BigFloat}, w::Complex{BigFloat}) = (operate!(-, z.re, w.re); operate!(-, z.im, w.im); z)
MA.operate!(::typeof(-), z::Complex{BigFloat}) = (operate!(-, z.re); operate!(-, z.im); z)

# multiplication

MA.promote_operation(::typeof(*), ::Type{Complex{BigFloat}}, ::Type{Complex{BigFloat}}) = Complex{BigFloat}

function MA.operate_to!(out::Complex{BigFloat}, ::typeof(*), a::Complex{BigFloat}, b::Complex{BigFloat})
    if (out === a === b)
        cbuf = cplx_buffer[][Threads.threadid()]
        operate_to!(cbuf, copy, a) # copy a into the buffer
        operate_to!(out.re, *, cbuf.re, cbuf.re)
        buffered_operate!(real_buffer[][Threads.threadid()], sub_mul, out.re, cbuf.im, cbuf.im)
        operate_to!(out.im, *, cbuf.re, cbuf.im)
        buffered_operate!(real_buffer[][Threads.threadid()], add_mul, out.im, cbuf.re, cbuf.im)
        return out
    elseif (out === a)
        cbuf = cplx_buffer[][Threads.threadid()]
        operate_to!(cbuf, copy, a) # copy a into the buffer
        a = copy(cbuf) # a is now the copied buffer -> de-aliased
    elseif (out === b)
        cbuf = cplx_buffer[][Threads.threadid()]
        operate_to!(cbuf, copy, b) # copy a into the buffer
        b = copy(cbuf) # a is now the copied buffer -> de-aliased
    end
    operate_to!(out.re, *, a.re, b.re)
    buffered_operate!(real_buffer[][Threads.threadid()], sub_mul, out.re, a.im, b.im)
    operate_to!(out.im, *, a.re, b.im)
    buffered_operate!(real_buffer[][Threads.threadid()], add_mul, out.im, a.im, b.re)
    return out
end

MA.operate!(::typeof(*), a::Complex{BigFloat}, b::Complex{BigFloat}) = operate_to!(a, *, a, b)

MA.operate_to!(out::Complex{BigFloat}, ::typeof(*), a::Complex{BigFloat}) = operate_to!(out, copy, a)
MA.operate!(::typeof(*), x::Complex{BigFloat}) = x  # unary + no-op

# Base.abs

function MA.operate_to!(o::BigFloat, ::typeof(abs2), x::Complex{BigFloat})
    operate_to!(o, *, x.re, x.re)
    buffered_operate!(real_buffer[][Threads.threadid()], add_mul, o, x.im, x.im)
    return o
end

function MA.operate!(::typeof(Base.abs2), x::BigFloat)
    o = BigFloat()
    return operate_to!(o, abs2, x)
end

# Base.conj

function MA.operate_to!(o::Complex{BigFloat}, ::typeof(Base.conj), x::Complex{BigFloat})
    operate_to!(o, copy, x)
    operate!(-, o.im)
    return o
end

function MA.operate!(::typeof(Base.conj), x::Complex{BigFloat})
    operate!(-, x.im)
    return x
end

# division

MA.promote_operation(::typeof(/), ::Type{Complex{BigFloat}}, ::Type{Complex{BigFloat}}) = Complex{BigFloat}
MA.promote_operation(::typeof(/), ::Type{Complex{BigFloat}}, ::Type{BigFloat}) = Complex{BigFloat}
MA.promote_operation(::typeof(/), ::Type{BigFloat}, ::Type{BigFloat}) = BigFloat

function MA.operate_to!(out::BigFloat, ::typeof(/), a::BigFloat, b::BigFloat)
    ccall((:mpfr_div, Base.MPFR.libmpfr), Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        out, a, b, Base.MPFR.ROUNDING_MODE[])
    return out
end

function MA.operate_to!(out::Complex{BigFloat}, ::typeof(/), a::Complex{BigFloat}, b::BigFloat)
    operate_to!(out.re, /, a.re, b) 
    operate_to!(out.im, /, a.im, b) 
    return out
end

function MA.operate_to!(out::Complex{BigFloat}, ::typeof(/), a::Complex{BigFloat}, b::Complex{BigFloat})
    bbar = cplx_buffer2[][Threads.threadid()]
    operate_to!(bbar, conj, b)
    m = real_buffer2[][Threads.threadid()]
    operate_to!(m, abs2, b)
    operate_to!(out, *, a, bbar)
    operate!(/, out, m)
    return out
end

function MA.operate!(::typeof(/), a::Complex{BigFloat}, b::BigFloat)
    operate!(/, a.re, b)
    operate!(/, a.im, b)
    return a
end

function MA.operate!(::typeof(/), a::Complex{BigFloat}, b::Complex{BigFloat})
    operate_to!(a, /, a, b)
    return a
end

# add_mul
# MA.promote_operation(::typeof(add_mul), ::Complex{BigFloat}, ::Complex{BigFloat}, ::Complex{BigFloat}) = Complex{BigFloat}

# function MA.operate_to!(out::Complex{BigFloat}, ::typeof(-), a::Complex{BigFloat}, b::Complex{BigFloat})
#     operate_to!(out.re, -, a.re, b.re)
#     operate_to!(out.im, -, a.im, b.im)
#     return out
# end

# MA.operate!(::typeof(-), z::Complex{BigFloat}) = (operate!(-, z.re); operate!(-, z.im); z)
