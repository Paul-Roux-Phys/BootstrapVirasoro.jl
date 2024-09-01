using JuliVirBootstrap
using Test

@testset "JuliVirBootstrap Tests" begin

    @testset "CFTData.jl" begin
        include("cft_data_tests.jl")
    end

    @testset "Conformal blocks" begin
        include("conformal_blocks/conformal_blocks_tests.jl")
    end

end
