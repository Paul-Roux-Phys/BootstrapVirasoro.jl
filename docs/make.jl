push!(LOAD_PATH, joinpath("..", "src"))

using BootstrapVirasoro
using Documenter

DocMeta.setdocmeta!(
    BootstrapVirasoro,
    :DocTestSetup,
    :(using BootstrapVirasoro),
    recursive = true
)

makedocs(
    sitename = "BootstrapVirasoro.jl",
    repo = "https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl",
    modules = [BootstrapVirasoro],
    format = Documenter.HTML(
        repolink = "https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl",
        edit_link = :commit
    ),
    doctest = true,
    pages = [
        "Home" => "index.md",
        "installation.md",
        "cft_data.md",
        "conformal_blocks.md",
        "bootstrap_equations.md",
        "reference.md"
    ],
    checkdocs=:export
)

deploydocs(;
    repo   = "https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl"
)
