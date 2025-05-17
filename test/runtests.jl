using BootstrapVirasoro,
    Test,
    Documenter

DocMeta.setdocmeta!(
    BootstrapVirasoro,
    :DocTestSetup,
    :(using BootstrapVirasoro),
    recursive=true
)
Documenter.doctest(BootstrapVirasoro)

@testset "CFTData.jl" begin
    include("cft_data.jl")
end

@testset "Conformal blocks" begin
    @testset "Residues" begin
        include("conformal_blocks/residues.jl")
    end

    @testset "Four point blocks" begin
        include("conformal_blocks/four_point.jl")
    end

    @testset "One point blocks" begin
        include("conformal_blocks/one_point.jl")
    end
end

# @testset "Bootstrap equations" begin
#     include("bootstrap_equations.jl")
# end
