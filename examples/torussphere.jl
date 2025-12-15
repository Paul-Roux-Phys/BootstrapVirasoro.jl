Pkg.activate("examples"); # activate the parent environment
println("Number of threads: $(Threads.nthreads())")
using BootstrapVirasoro, BootstrapVirasoro.LoopModels, BenchmarkTools, Test

const Sig = Channels{Rational}

# Spectrum of O(n) model (sphere and torus central charges)
c = CC(β = sqrt(big"2") / (big"0.8" + big"0.1" * im))
c_S2 = CC(β = 1 / (big"0.8" + big"0.1" * im))
ndiag_indices = [
    (r, s) for r = (1//2):(1//2):20 for s = (-1+1//(2r)):(1//(2r)):1 if r * s % 1 // 2 == 0
]
ndiag_indices_S2 = [(r, s) for (r, s) in ndiag_indices if r * s % 1 == 0]
diag_field = Field(c, r = 0, s = big"0.4" + big"0.1" * im)
diag_field_S2 = Field(c, r = 0, s = diag_field.s)
@test diag_field.dims.left.P ≈ diag_field_S2.dims.left.P
fields = vcat([Field(c, r = r, s = s) for (r, s) in ndiag_indices], diag_field)
fields_S2 =
    vcat([Field(c_S2, r = r, s = s) for (r, s) in ndiag_indices_S2], diag_field_S2)

# Determine the parity of the number of legs in 4pt channels
function chan_parities(co::Correlation4)
    V1, V2, V3, V4 = co.fields
    chan_parities =
        Channels{Rational}((V1.r + V2.r) % 1, (V1.r + V4.r) % 1, (V1.r + V3.r) % 1)
end

# Compute series of blocks, keeping only one of each (field, reflected field) pair
# given indices of the external fields.
# torus version
function precompute_blocks(
    index::NTuple{2,<:Number},
    fields;
    parity = 0,
    precision = 10,
)
    setprecision(BigFloat, floor(Int, 1.2 * precision), base = 10)
    Δmax = floor(Int, 1.5 * precision)
    co = Correlation(Field(c, r = index[1], s = index[2]), Δmax)
    parity != 0 && (fields = filter(V -> V.diagonal || V.s >= 0, fields))
    f = V -> IBlock(co, V, parity = parity)
    spec = ChanSpec(co, fields, f)
    return @channels spec
end
# sphere version
function precompute_blocks(indices, fields; parity = 0, precision = 10)
    setprecision(BigFloat, floor(Int, 1.2 * precision), base = 10)
    Δmax = floor(Int, 1.5 * precision)
    co = Correlation([Field(c_S2, r = r, s = s) for (r, s) in indices], Δmax)
    parity != 0 && (fields = filter(V -> V.diagonal || V.s >= 0, fields))
    fs = Channels(
        chan -> (V -> IBlock(co, chan, V, parity = parity)),
    )
    Vs = @channels filter(V -> V.r % 1 == chan_parities(co)[chan], fields)
    return @channels ChannelSpectrum(co, chan, Vs[chan], fs[chan])
end

# prefactors for torus-sphere correspondence
function prefA(τ, co_s, chan, lr)
    x_cache = PosCache(xfromq(exp(im * (π * τ))), co_s[lr], chan)
    return inv(etaDedekind(τ) * x_cache.prefactor)
end

ncprefA(τ, co_s, chan) =
    prefA(τ, co_s, chan, :left) * prefA(-conj(τ), co_s, chan, :right)

# Solve crossing symmetry for given signature, show solution
function solve(specs, signature)
    specs = @channels filter(V -> V.r >= signature[chan], specs[chan])
    system = solve_bootstrap(specs)
    show(stdout, system.str_cst, rmax = 3)
    return system
end

relative_error(a, b) = abs((a - b) / (a + b))

setprecision(BigFloat, 40, base = 10)
c = CentralCharge(β = big"1.2" + big"0.1" * im)
c_S2 = CentralCharge(β = c.β / sqrt(big(2)))
fields = vcat([Field(c, r = r, s = s) for (r, s) in ndiag_indices], diag_field)
diag_field = Field(c, r = 0, s = (big"0.4" + big"0.1" * im) / 2)
diag_field_S2 = Field(c, r = 0, s = diag_field.s)
fields_S2 =
    vcat([Field(c_S2, r = r, s = s) for (r, s) in ndiag_indices], diag_field_S2)
Δmax = 40

ind = (1, 0)
V1 = Field(c, r = ind[1], s = ind[2])
co = Correlation(V1, Δmax)

ind_S2 = [(0, 1 // 2), (ind[1], ind[2] / 2), (0, 1 // 2), (0, 1 // 2)]
co_S2 = Correlation([Field(c_S2, r = r, s = s) for (r, s) in ind_S2], Δmax)

blocks = precompute_blocks(ind, fields, precision = 5, parity = 1)
blocks_S2 = precompute_blocks(ind_S2, fields_S2, precision = 8, parity = 1)

sys_S2 = BootstrapSystem(blocks_S2)
sys = BootstrapSystem(blocks, pos = τfromx.(sys_S2.positions))

# for V in keys(blocks.s.blocks)
#     V_S2 = Field(c_S2, r = 2 * V.r, s = V.s)
#     if V_S2 in keys(blocks_S2.s.blocks)
V = Field(c, r=3//2, s=1//3)
V_S2 = Field(c_S2, r=2*V.r, s=V.s)
b = blocks.s.blocks[V]
b_S2 = @channels blocks_S2[chan].blocks[V_S2]

P = b.chan_field.dims.left.P
Pbar = b.chan_field.dims.right.P
τ = big"0.36" + big"1.14" * im
x = xfromq(exp(im * (π * τ)))

# for chan in BootstrapVirasoro.CHANNELS
chan = :u
    t = BootstrapVirasoro.modular_param(chan, τ)
    X = BootstrapVirasoro.crossratio(chan, x)
    err = relative_error(
        b(t) * 16^(2P^2 + 2Pbar^2),
        b_S2[chan](X) * ncprefA(t, co_S2, chan),
    )

    if err > 1e-10
        println("$V, $chan, $err")
    end
# end
#     end
# end

evaluate_blocks!(sys)
evaluate_blocks!(sys_S2)
compute_linear_system!(sys)
compute_linear_system!(sys_S2)

solve!(sys)
solve!(sys_S2)

sort(
    [
        (V, ncprefA(sys.positions[1], co_S2, :s) * b(sys_S2.positions[1])) for
        (V, b) in sys_S2.spectra.s.blocks
    ],
    by = b -> real(b[2]),
)
sort(
    [
        (V, 16^(2V.dims.left.P^2 + 2V.dims.right.P^2) * b(sys.positions[1])) for
        (V, b) in sys.spectra.s.blocks
    ],
    by = b -> real(b[2]),
)

solve(blocks, Sig(1 // 2, 0, 1 // 2))

V = Field(c, r = 3 // 2, s = 2 // 3)
V_S2 = Field(c_S2, r = 3, s = 2 // 3)
Block(co, Field(c, r = 3 // 2, s = 2 // 3), interchiral = true, parity = 1)(
    sys.positions[1],
) * 16^(2V.dims.left.P^2 + 2V.dims.right.P^2)
Block(co_S2, :s, Field(c_S2, r = 3, s = 2 // 3), interchiral = true, parity = 1)(
    sys_S2.positions[1],
) * ncprefA(sys.positions[1], co_S2, :s)
