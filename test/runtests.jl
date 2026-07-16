using BootstrapVirasoro, Test

@testset "CFTData.jl" begin
    include("cft_data.jl")
end

@testset "Four point blocks" begin
    include("four_point.jl")
end

@testset "One point blocks" begin
    include("one_point.jl")
end

@testset "Loop models submodule" begin
    include("LoopModels.jl")
end

@testset "Examples" begin
    include("../examples/compute_blocks.jl")
    include("../examples/bootstrap.jl")

    for chan in (:s, :t)
        for (V, (val, err)) in sol[chan].dict
            @test abs(err) < 1e-20
        end
    end
end
