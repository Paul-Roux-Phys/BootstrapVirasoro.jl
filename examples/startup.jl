using Pkg: Pkg
let
    pkgs = [
        "Revise", "IJulia", "LanguageServer",
        "JuliaFormatter", "Plots", "OhMyREPL",
        "BenchmarkTools"
    ]
    for pkg in pkgs
        if Base.find_package(pkg) === nothing
            Pkg.add(pkg)
        end
    end
end

try
    using Revise
    ENV["JULIA_REVISE"] = "auto"
catch e
    @warn "Error initializing Revise" exception = (e, catch_backtrace())
end

try
    using BenchmarkTools
catch e
    @warn "Error initialising BenchmarkTools" exception = (e, catch_backtrace())
end

atreplinit() do repl
    try
        @eval using OhMyREPL
    catch e
        @warn "Error initializing OhMyRepl" exception = (e, catch_backtrace())
    end
end

if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
