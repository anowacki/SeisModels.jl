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

        @testset "ek137" begin
            @test EK137 isa LinearLayeredModel
            @test surface_radius(EK137) == 6371.0
            @test length(EK137.r) == 138
            @test vp(EK137, 0) ≈ 11.2622
            @test !isanisotropic(EK137)
            @test hasdensity(EK137)
            @test !hasattenuation(EK137)
            @test_throws ArgumentError Qμ(EK137, 1000)
            @test density(EK137, 0, depth=true) ≈ 2.6510
        end

        @testset "iasp91" begin
            @test surface_radius(IASP91) == 6371.0
            @test length(IASP91.r) == 11
            @test !isanisotropic(IASP91)
            @test !hasdensity(IASP91)
            @test !hasattenuation(IASP91)
            @test_throws ArgumentError mass(IASP91)
        end

        @testset "PREM" begin
            @test surface_radius(PREM) == 6371.0
            @test length(PREM.r) == 13
            @test isanisotropic(PREM)
            @test vp(PREM, 0) ≈ 11.2622
            @test hasattenuation(PREM)
            @test hasreffrequency(PREM)
            @test reffrequency(PREM) == 1.0
            @test Qμ(PREM, 1000) ≈ 84.6
            @test eta(PREM, 1000) == 1.0
            @test moment_of_inertia(PREM)/(
                mass(PREM)*surface_radius(PREM)^2*1e6) ≈ 0.3308 atol=0.0001
            @test gravity(PREM, surface_radius(PREM)) ≈ 9.81 atol=0.02
        end

        @testset "PREM_NOOCEAN" begin
            @test surface_radius(PREM_NOOCEAN) == surface_radius(PREM)
            @test length(PREM_NOOCEAN.r) == length(PREM.r) - 1
            @test isanisotropic(PREM_NOOCEAN)
            @test vp(PREM_NOOCEAN, 0) ≈ 11.2622
            @test hasattenuation(PREM_NOOCEAN)
            @test hasreffrequency(PREM_NOOCEAN)
            @test reffrequency(PREM_NOOCEAN) == 1.0
            @test Qμ(PREM_NOOCEAN, 0, depth=true) == 600.0
            @test eta(PREM_NOOCEAN, 0, depth=true) == 1.0
            @test vp(PREM_NOOCEAN, 0, depth=true) == 5.8
            @test vsv(PREM_NOOCEAN, 0, depth=true) == 3.2
        end

        @testset "SC_REM" begin
            sc_rem_models = (
                SC_REM_25_LOW, SC_REM_25_HIGH, SC_REM_50_LOW, SC_REM_50_HIGH,
                SC_REM_75_LOW, SC_REM_75_HIGH
            )
            for m in sc_rem_models
                @test surface_radius(m) == 6371.0
                @test !isanisotropic(m)
                @test hasdensity(m)
                @test hasattenuation(m)
                I_MR2 = moment_of_inertia(m)/(mass(m)*surface_radius(m)^2*1e6)
                # Eyeballed from Figure 4 of Kemper et al. (2023)
                @test 0.3297 <= I_MR2 <= 0.331
            end
        end

        @testset "STW105" begin
            @test surface_radius(STW105) == 6371.0
            @test isanisotropic(STW105)
            @test hasattenuation(STW105)
            @test hasdensity(STW105)
            @test length(STW105.r) == 750
            @test mass(STW105) ≈ 5.974e24 atol=0.0005e24
            @test moment_of_inertia(STW105)/(
                mass(STW105)*surface_radius(STW105)^2*1e6) ≈ 0.3308 atol=0.0001
            @test eta(STW105, 100, depth=true) ≈ 0.85 atol=0.1
            @test density(STW105, 150.127) ≈ 13.08357 atol=0.00001
            @test vsh(STW105, 6291.0001) ≈ 4.57246
            @test Qμ(STW105, 3482) == 355.0
            @test Qκ(STW105, 0) == 1327.6
        end
    end

    @testset "Moon models" begin
        @testset "Weber et al. (Science, 2011)" begin
            @test surface_radius(MOON_WEBER_2011) == 1737.1
            @test !isanisotropic(MOON_WEBER_2011)
            @test hasdensity(MOON_WEBER_2011)
            @test !hasattenuation(MOON_WEBER_2011)
            @test mass(MOON_WEBER_2011) ≈ 7.35e22 atol=0.01e22
            @test moment_of_inertia(MOON_WEBER_2011) / 
                (mass(MOON_WEBER_2011) *
                    surface_radius(MOON_WEBER_2011)^2 * 1e6) ≈ 0.394 atol=0.01
        end
    end
end
