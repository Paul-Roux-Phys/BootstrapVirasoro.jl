@testset "Residues" begin
    include("zamolodchikov_residues_tests.jl")
end

@testset "Four point blocks" begin
    include("four_point_blocks_tests.jl")
end

@testset "One point blocks" begin
    include("one_point_blocks_tests.jl")
end
