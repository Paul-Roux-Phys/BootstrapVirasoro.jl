using JuliVirBootstrap
using Documenter

DocMeta.setdocmeta!(JuliVirBootstrap, :DocTestSetup, :(using JuliVirBootstrap); recursive=true)

makedocs(;
    modules=[JuliVirBootstrap],
    authors="Paul_Roux",
    repo="https://gitlab.com/Paul Roux/JuliVirBootstrap.jl/blob/{commit}{path}#{line}",
    sitename="JuliVirBootstrap.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Paul Roux.gitlab.io/JuliVirBootstrap.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
