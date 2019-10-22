using Test
using EarthModels

@testset "1D models" begin
    @testset "Find layer" begin
        @test_throws DomainError EarthModels.findlayer(PREM, -1)
        @test_throws DomainError EarthModels.findlayer(PREM, 1e6)
    end
end