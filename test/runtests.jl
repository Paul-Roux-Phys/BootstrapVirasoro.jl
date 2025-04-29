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
        include("cft_data_tests.jl")
    end

    @testset "Conformal blocks" begin
        @testset "Residues" begin
            include("conformal_blocks/zamolodchikov_residues_tests.jl")
        end

        @testset "Four point blocks" begin
            include("conformal_blocks/four_point_blocks_tests.jl")
        end

        @testset "One point blocks" begin
            include("conformal_blocks/one_point_blocks_tests.jl")
        end
    end

    # @testset "Bootstrap equations" begin
    #     include("bootstrap_equations.jl")
    # end

end
