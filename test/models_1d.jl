using Test
using SeisModels

@testset "1D models" begin
    @testset "Broadcasting" begin
        @test vp.((PREM, AK135), 1000) == (11.10541122721928, 11.103425443786982)
        @test vs.(PREM, 0:3) ≈ [3.6678, 3.667799890427708, 3.6677995617108317, 3.6677990138493715]
    end

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
                # Evaluation via property symbol
                for (f, sym) in zip((vp, vs, density, vph, vpv, vsh, vsv, eta, Qμ,  Qκ),
                    (:vp, :vs, :density, :vph, :vpv, :vsh, :vsv, :eta, :Qμ,  :Qκ))
                    @test f(m, 0.9) == evaluate(m, sym, 0.9)
                end
            end
        end

        # TODO: Add tests for LinearLayeredModel and PREMPolyModel
    end

    # Test for equivalence of depth versions of functions
    @testset "Depth" begin
        let n = 5, radii = sort(rand(n)), R = radii[end],
                (Vp, Vs, rho, Vpv, Vph, Vsv, Vsh, Eta, Qkappa, Qmu) = rand.(repeat([n], 10))
            for model_type in (SteppedLayeredModel, LinearLayeredModel, PREMPolyModel)
                m = if model_type <: PREMPolyModel
                    model_type(R, n, radii, Vp', Vs', rho', true, Vpv', Vph', Vsv', Vsh', Eta',
                        true, Qkappa', Qmu', NaN)
                else
                    model_type(R, n, radii, Vp, Vs, rho, true, Vpv, Vph, Vsv, Vsh, Eta,
                        true, Qkappa, Qmu)
                end
                r = R*0.99*rand()
                z = R - r
                for func in (vp, vs, density, vpv, vph, vsv, vsh, eta, Qκ, Qμ,
                        pressure, bulk_modulus, shear_modulus,
                        youngs_modulus, poissons_ratio, mass, surface_mass)
                    @test func(m, z, depth=true) ≈ func(m, r)
                end
                @test evaluate(m, :vp, z, depth=true) ≈ evaluate(m, :vp, r)
            end
        end
    end
end
