using BootstrapVirasoro,
    Test,
    Documenter

@testset "BootstrapVirasoro Tests" begin

    # DocMeta.setdocmeta!(
    #     BootstrapVirasoro,
    #     :DocTestSetup,
    #     :(using BootstrapVirasoro),
    #     recursive=true
    # )
    # Documenter.doctest(BootstrapVirasoro)

    @testset "CFTData.jl" begin
        include("cft_data/cft_data_tests.jl")
    end

    @testset "Conformal blocks" begin
        include("conformal_blocks/conformal_blocks_tests.jl")
    end

end
