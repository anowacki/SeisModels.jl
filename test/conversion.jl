using Test, SeisModels

@testset "Conversion" begin
    # Default spacing
    let m = LinearLayeredModel(PREM)
        @test m isa LinearLayeredModel
        @test maximum(diff(m.r)) == 20
        @test minimum(diff(m.r)) == 0
        @test vp(m, 1000) ≈ vp(PREM, 1000) rtol=0.0001
    end
    # Specified spacing
    let m = LinearLayeredModel(IASP91, 30)
        @test maximum(diff(m.r)) == 30
        @test vs(m, 3600) ≈ vs(IASP91, 3600) rtol=0.0001
    end
end