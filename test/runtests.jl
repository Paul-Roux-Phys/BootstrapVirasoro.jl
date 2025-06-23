using BootstrapVirasoro, Test, Documenter

DocMeta.setdocmeta!(
    BootstrapVirasoro,
    :DocTestSetup,
    :(using BootstrapVirasoro),
    recursive = true,
)
Documenter.doctest(BootstrapVirasoro)

@testset "CFTData.jl" begin
    include("cft_data.jl")
end

@testset "Block residues" begin
    include("residues.jl")
end

@testset "Four point blocks" begin
    include("four_point.jl")
end

@testset "One point blocks" begin
    include("one_point.jl")
end

# @testset "Bootstrap equations" begin
#     include("bootstrap_equations.jl")
# end
