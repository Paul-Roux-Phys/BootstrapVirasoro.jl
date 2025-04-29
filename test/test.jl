n = 5
setprecision(BigFloat, 256)
a = rand(Complex{BigFloat}, n)
qs = rand(Complex{BigFloat}, 10)

buf = Complex{BigFloat}(0)
buf_mul = Complex{BigFloat}(0)

[evalpoly(q, a) for q in qs]
[evalpoly_buf(q, a) for q in qs]
[evalpoly_buf!(buf_mul, buf, q, a) for q in qs]
evalpoly_buf(qs, a)

@btime [evalpoly(q, a) for q in qs]
@btime [evalpoly_buf(q, a) for q in qs]
@btime [evalpoly_buf!(buf_mul, buf, q, a) for q in qs]
@btime evalpoly_buf(qs, a) # threaded implementation

println(@allocated [evalpoly(q, a) for q in qs])
println(@allocated [evalpoly_buf(q, a) for q in qs])
println(@allocated [evalpoly_buf!(buf_mul, buf, q, a) for q in qs])
println(@allocated evalpoly_buf(qs, a))

