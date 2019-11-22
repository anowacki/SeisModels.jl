using SeisModels, Test

@testset "IO" begin
    mktempdir() do dir
        file = joinpath(dir, "model")
        write_mineos(AK135, file)
        ak135 = read_mineos(file)
        @test AK135 ≈ ak135
    end
end
