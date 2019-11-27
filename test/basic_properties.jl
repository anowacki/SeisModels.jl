using Test
using SeisModels

@testset "Basic properties" begin
    let m = SteppedLayeredModel(1, 2, [0.5, 1], [1, 2], [3, 4], [5, 6], false,
            [], [], [], [], [], true, [10, 20], [100, 200])
        @test !isanisotropic(m)
        @test surface_radius(m) == 1.0
        @test radius(m, 0.25) == 0.75
        @test depth(m, 0.25) == 0.75
        @test hasattenuation(m)
        @test hasdensity(m)
    end
    let a = 1, n = 2, arr = rand(2), m = SteppedLayeredModel(a, n, [0.005 + 0.99rand(), 1],
            arr, arr, [], true, arr, arr, arr, arr, arr, false, [], [])
        @test isanisotropic(m)
        @test !hasattenuation(m)
        @test !hasdensity(m)
    end
    @test hasreffrequency(PREM)
    @test reffrequency(PREM) == 1.0
    let m = PREMPolyModel([getfield(PREM, f)
                           for f in fieldnames(PREMPolyModel)][1:end-1]..., NaN)
        @test !hasreffrequency(m)
        @test_throws ArgumentError reffrequency(m)
    end
end

@testset "Model comparison" begin
    @test IASP91 == deepcopy(IASP91)
    @test PREM != AK135
    let prem′ = deepcopy(PREM)
        @test prem′ ≈ PREM
        prem′.vp .= 1.00000001 .*PREM.vp
        @test prem′ ≈ PREM
        prem′.vp .= 1.000001 .* PREM.vp
        @test !(prem′ ≈ PREM)
    end
end
