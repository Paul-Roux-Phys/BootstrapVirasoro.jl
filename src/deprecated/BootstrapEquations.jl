using Pkg; Pkg.activate(".")
using JuliVirBootstrap

help(Field)

function evaluate_block(positions, Nmax, corr, block)
    res = zeros(length(positions))
    threads.@Threads for (i,pos) in enumerate(positions)
        res[i] = G(corr, block, pos)
    end
end

