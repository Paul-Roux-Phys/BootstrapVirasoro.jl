Pkg.activate(".") # activate the parent environment
using BootstrapVirasoro, BenchmarkTools
n_threads = Threads.nthreads()
if n_threads == 1
    println("You are using a single thread. Consider starting julia with more threads, for instance by setting
    the environment variable `export JULIA_NUM_THREADS=auto`")
end

setprecision(BigFloat, 35, base=10)
c = CC(β=big"0.8" + big"0.1"*im)
Δmax = 40.

function LoopFields(model)
    model === :On && return vcat(
        [Field(c, r=r, s=s) for r in 1:15 for s in -1+1//r:1//r:1],
        [Field(c, r=r, s=s) for r in 1//2:1:15 for s in -1+1//(2r):1//r:1 if (r*s)%1 == 0],
        [Field(c, r=1, s=s, diagonal=true) for s in 1:2:15]
    )
    model === :PSUn && return vcat(
        [Field(c, r=r, s=s) for r in 1:15 for s in -1+1//r:1//r:1],
        [Field(c, r=1, s=s, diagonal=true) for s in 1:15]
    )
    model === :Potts && return vcat(
        [Field(c, r=r, s=s) for r in 2:15 for s in -1+1//r:1//r:1],
        [Field(c, r=0, s=s, diagonal=true) for s in 1//2:1:3//2],
        [Field(c, r=1, s=s, diagonal=true) for s in 1:15]
    )
end

LoopSpectrum(model, Δmax) = Spectrum(LoopFields(model), Δmax, interchiral=true);

indices = ((1//2, 0), (1//2, 0), (1, 0), (1, 0))
fields = [Field(c, r=r, s=s) for (r, s) in indices]
co = Correlation(fields..., Δmax=Δmax)
println("time to compute residues:")
Δmax = 6.
@btime Correlation(fields..., Δmax=Δmax)
SPSUn = LoopSpectrum(:On, Δmax)
println("time to compute the spectrum:")
@btime LoopSpectrum(:On, Δmax)
Schan = ChannelSpectra(co, SPSUn)
println("time to precompute the block coefficients:")
@btime ChannelSpectra(co, SPSUn)
PrecomputedSystem = BootstrapSystem(Schan)
println("time to compute all blocks evaluated at all positions:")
@btime BootstrapSystem(Schan)

function BootstrapSystem(b; signature=(s=0, t=0, u=0), diags=(s=nothing, t=nothing, u=nothing))
    new = deepcopy(b)
    for chan in (:s, :t, :u), V in b.spectra[chan].fields
        if (isdiagonal(V) && signature[chan] > 0) || V.r < signature[chan]
            remove!(new, chan, V)
        end
        if !isnothing(diags[chan])
            for V in diags[chan]
                add!(new, chan, V)
            end
        end
    end
    return new
end

println("time to form the system for a given signature:")
sys = BootstrapSystem(PrecomputedSystem, signature=(s=1, t=1 // 2, u=3 // 2))
@btime BootstrapSystem(PrecomputedSystem, signature=(s=1, t=1 // 2, u=3 // 2))
println("time to form the matrix of the system")
compute_linear_system!(sys)
@btime compute_linear_system!(sys)
println("time to solve the system")
solve!(sys)
@btime solve!(sys)