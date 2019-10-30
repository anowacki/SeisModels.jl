using Test
using SeisModels

@testset "1D models" begin
    @testset "Find layer" begin
        @test_throws DomainError SeisModels.findlayer(PREM, -1)
        @test_throws DomainError SeisModels.findlayer(PREM, 1e6)
    end

    @testset "Evaluation" begin
        @testset "SteppedLayeredModel" begin
            let m = SteppedLayeredModel(1, 2, [0.25, 1], [2, 1], [0.5, 1], [4, 3],
                    true, [1, 2], [3, 4], [5, 6], [7, 8], [9, 10],
                    true, [100, 200], [300, 400])
                @test vp(m, 0) == 2
                @test vp(m, 0.26) == 1
                @test vs(m, 0.1) == 0.5
                @test vs(m, 1) == 1
                @test density(m, 0.1) == 4
                @test density(m, 0.9) == 3
                @test vph(m, 0.2) == 1
                @test vph(m, 0.25) == 2
                @test vpv(m, 0) == 3
                @test vpv(m, 1) == 4
                @test vsh(m, 0.24) == 5
                @test vsh(m, nextfloat(0.25)) == 6
                @test vsv(m, prevfloat(0.25)) == 7
                @test vsv(m, nextfloat(0.25)) == 8
                @test eta(m, prevfloat(0.25)) == 9
                @test eta(m, nextfloat(0.25)) == 10
                @test Qmu(m, 0) == Qμ(m, 0) == 100
                @test Qμ(m, 1) == 200
                @test Qκ(m, 0.2) == Qkappa(m, 0.2) == 300
                @test Qκ(m, 0.3) == 400
                @test mass(m, surface_radius(m)) == 4/3*π*(0.25e3^3*4e3 + (1e3^3 - 0.25e3^3)*3e3)
                @test gravity(m, 0.75) ≈ 6.67428e-11*4/3*π*(0.25e3^3*4e3 + (0.75e3^3 - 0.25e3^3)*3e3)/0.75e3^2
                @test shear_modulus(m, 1) == 3e3*1e3^2/1e9
                @test bulk_modulus(m, 0.5) == (3e3*1e3^2/1e9 - 4/3*shear_modulus(m, 0.5))
                @test poissons_ratio(m, 0) == 0.4666666666666667
            end
        end
    end
end
