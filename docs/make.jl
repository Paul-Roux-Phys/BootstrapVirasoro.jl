push!(LOAD_PATH, joinpath("..", "src"))

using Documenter
using Documenter.Remotes: GitHub
using BootstrapVirasoro

DocMeta.setdocmeta!(
        BootstrapVirasoro,
        :DocTestSetup,
        :(using BootstrapVirasoro),
)

makedocs(
        sitename="BootstrapVirasoro.jl",
        repo=GitHub("Paul-Roux-Phys", "BootstrapVirasoro.jl"),
        modules=[BootstrapVirasoro],
        format=Documenter.HTML(
                repolink="https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl",
                edit_link=:commit,
        ),
        pages=[
                "Home" => "index.md",
                "Installation" => "installation.md",
                "Conformal blocks" => "conformalblocks.md",
                "Bootstrap equations" => "bootstrapeqs.md",
        ],
        checkdocs=:export,
)

deploydocs(
        repo=GitHub("Paul-Roux-Phys", "BootstrapVirasoro.jl"),
        devbranch="main",
)
