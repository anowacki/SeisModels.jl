using Test
using SeisModels

@testset "Inbuilt models" begin
    @testset "Earth models" begin
        @testset "ak135" begin
            @test surface_radius(AK135) == 6371.0
            @test length(AK135.r) == 136
            @test vp(AK135, 0) ≈ 11.2622
            @test !isanisotropic(AK135)
            @test !hasattenuation(AK135)
            @test_throws ArgumentError Qμ(AK135, 1000)
        end

        @testset "iasp91" begin
            @test surface_radius(IASP91) == 6371.0
            @test length(IASP91.r) == 11
            @test !isanisotropic(IASP91)
            @test !hasdensity(IASP91)
            @test !hasattenuation(IASP91)
        end

        @testset "PREM" begin
            @test surface_radius(PREM) == 6371.0
            @test length(PREM.r) == 13
            @test isanisotropic(PREM)
            @test vp(PREM, 0) ≈ 11.2622
            @test hasattenuation(PREM)
            @test Qμ(PREM, 1000) ≈ 84.6
            @test moment_of_inertia(PREM)/(
                mass(PREM, surface_radius(PREM))*surface_radius(PREM)^2*1e6) ≈ 0.3308 atol=0.0001
            @test gravity(PREM, surface_radius(PREM)) ≈ 9.81 atol=0.02
        end
    end

    @testset "Moon models" begin
        @testset "Weber et al. (Science, 2011)" begin
            @test surface_radius(MOON_WEBER_2011) == 1737.1
            @test !isanisotropic(MOON_WEBER_2011)
            @test hasdensity(MOON_WEBER_2011)
            @test !hasattenuation(MOON_WEBER_2011)
            @test mass(MOON_WEBER_2011, surface_radius(MOON_WEBER_2011)) ≈ 7.35e22 atol=0.01e22
            @test moment_of_inertia(MOON_WEBER_2011) / 
                (mass(MOON_WEBER_2011, surface_radius(MOON_WEBER_2011)) *
                    surface_radius(MOON_WEBER_2011)^2 * 1e6) ≈ 0.394 atol=0.01
        end
    end
end
