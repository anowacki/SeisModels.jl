using Test
using SeisModels

@testset "1D models" begin
    @testset "Find layer" begin
        @test_throws DomainError SeisModels.findlayer(PREM, -1)
        @test_throws DomainError SeisModels.findlayer(PREM, 1e6)
    end
end
