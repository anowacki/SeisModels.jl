using Test, SeisModels

@testset "Conversion" begin
    @testset "To LinearLayeredModel" begin
        # PREMPolyModel
        let m = LinearLayeredModel(PREM)
            @test m isa LinearLayeredModel
            # Default spacing
            @test maximum(diff(m.r)) == 20
            @test minimum(diff(m.r)) == 0
            # Specified spacing
            @test maximum(diff(LinearLayeredModel(IASP91, 30).r)) == 30
            @test vp(m, 1000) ≈ vp(PREM, 1000) rtol=0.0001
            @test vs(m, 3600) ≈ vs(PREM, 3600) rtol=0.0001
        end
        # SteppedLayeredModel
        let m = MOON_WEBER_2011, m′ = LinearLayeredModel(m)
            @test m′ isa LinearLayeredModel
            @test vp(m, 100) ≈ vp(m′, 100) rtol=0.0001
        end
    end

    @testset "To PREMPolyModel" begin
        # SteppedLayeredModel
        @test_throws ArgumentError PREMPolyModel(MOON_WEBER_2011, -1)
        @test reffrequency(PREMPolyModel(MOON_WEBER_2011, fref=10.0)) == 10.0
        let m = MOON_WEBER_2011, m′ = PREMPolyModel(MOON_WEBER_2011)
            for property in (:vp, :vs, :density)
                @test all(r -> evaluate(m, property, r) ≈ evaluate(m′, property, r),
                    range(0, stop=surface_radius(m), length=10_000))
            end
            @test mass(m) ≈ mass(m′)
            @test pressure(m, 0) ≈ pressure(m′, 0)
        end
        # LinearLayeredModel
        @test_throws ArgumentError PREMPolyModel(AK135, 0)
        @test_throws ArgumentError PREMPolyModel(AK135, -1)
        @test reffrequency(PREMPolyModel(AK135, fref=0.1)) == 0.1
        let m = AK135, m′ = PREMPolyModel(AK135)
            for property in (:vp, :vs, :density)
                @test all(r -> evaluate(m, property, r) ≈ evaluate(m′, property, r),
                    range(0, stop=surface_radius(m), length=10_000))
            end
            @test mass(m) ≈ mass(m′)
            @test pressure(m, 0) ≈ pressure(m′, 0)
        end
    end

    @testset "To SteppedLayeredModel" begin
        @test_throws ArgumentError SteppedLayeredModel(IASP91, 0)
        for m in (PREM, AK135)
            m′ = SteppedLayeredModel(m)
            @test m′ isa SteppedLayeredModel
            radii = range(0, stop=surface_radius(m), length=10_000)
            for property in (:vp, :vs, :density)
                @test all(r -> ≈(evaluate.((m, m′), property, r)..., rtol=0.01),
                    radii)
            end
            @test mass(m) ≈ mass(m′) rtol=0.001
            @test pressure(m, 0) ≈ pressure(m′, 0) rtol=0.001
        end
        @test maximum(diff(SteppedLayeredModel(PREM, 20).r)) == 20
    end

    @testset "Roundtrip" begin
        let rtol = 0.01
            for T in (PREMPolyModel, LinearLayeredModel, SteppedLayeredModel)
                for m in (PREM, AK135, MOON_WEBER_2011)
                    m isa T && continue # Skip identity conversions
                    m′ = T(m)
                    m″ = typeof(m)(m′)
                    models = (m, m″)
                    radii = range(0, stop=surface_radius(m), length=10_000)
                    for property in (:vp, :vs, :density)
                        @test all(r -> ≈(evaluate.(models, property, r)..., rtol=rtol),
                            radii)
                    end
                    @test ≈(mass.(models)..., rtol=rtol)
                    @test ≈(gravity.(models, surface_radius(m))..., rtol=rtol)
                    @test ≈(pressure.(models, 0)..., rtol=rtol)
                end
            end
        end
    end
end
