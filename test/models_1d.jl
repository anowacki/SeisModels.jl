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

        @testset "PREMPolyModels" begin
            # Velocities default to reference frequency
            for f in (vp, vs, vpv, vsv, vph, vsh)
                @test f(PREM, 5000) == f(PREM, 5000, freq=reffrequency(PREM))
            end
            # Value from separate implementation at https://github.com/andreww/prem4derg
            @test vs(PREM, 1000, freq=0.1) == 3.5274008549889566
        end

        # TODO: Add tests for LinearLayeredModel and PREMPolyModel
    end

    # Keyword constructors
    @testset "Construction" begin
        let n = 5, r = [0; sort(10_000rand(n-1))], Vp = rand(n), Vs = 0.5rand(n),
                Rho = rand(n), Vpv = 1.01.*Vp, Vph = 0.99.*Vp,
                Vsv = 0.99.*Vs, Vsh = 1.01.*Vs, Eta = 1 .+ 0.1rand(n),
                qμ = rand(n), qκ = rand(n), R = last(r), fref = rand(),
                null = Array{Float64}(undef, 0, 0)
            for model_type in (LinearLayeredModel, SteppedLayeredModel)
                # Arrays wrong length
                @test_throws ArgumentError model_type(r=r, vp=[Vp; 1])
                # r must be provided
                @test_throws UndefKeywordError model_type(vp=Vp, vs=Vs)
                # r must increase
                @test_throws ArgumentError model_type(r=[3, 2, 1])
                # Equivalence of default and keyword constructors
                @test model_type(R, n, r, Vp, Vs, Rho, true, Vph, Vpv,
                    Vsh, Vsv, Eta, true, qμ, qκ) ==
                    model_type(r=r, vp=Vp, vs=Vs, density=Rho, vph=Vph,
                        vpv=Vpv, vsh=Vsh, vsv=Vsv, Qμ=qμ, Qκ=qκ, eta=Eta)
                @test model_type(R, n, r, Vp, Vs, [], false, [],
                                         [], [], [], [], true, qμ, qκ) ==
                    model_type(r=r, vp=Vp, vs=Vs, Qμ=qμ, Qκ=qκ)
                @test model_type(R, n, r, Vp, Vs, Rho, true,
                                         Vph, Vpv, Vsh, Vsv, Eta, false, [], []) ==
                    model_type(r=r, vp=Vp, vs=Vs, density=Rho,
                                       vph=Vph, vpv=Vpv, vsh=Vsh, vsv=Vsv, eta=Eta)
                # ASCII arguments are the same and take precedence
                @test model_type(r=r, Qmu=qμ, Qkappa=qκ) ==
                    model_type(r=r, Qμ=qμ, Qκ=qκ) ==
                    model_type(r=r, Qμ=Vs, Qκ=Vp, Qmu=qμ, Qkappa=qκ)
            end
            ## Model-specific tests
            # LinearLayeredModel
            # r must start at 0
            @test_throws ArgumentError LinearLayeredModel(r=r[2:end])

            # PREMPolyModel
            @test_throws ArgumentError PREMPolyModel(r=r, vp=[Vp; 1])
            @test_throws UndefKeywordError PREMPolyModel(vp=Vp, vs=Vs)
            # r must increase
            @test_throws ArgumentError PREMPolyModel(r=[3, 2, 1])
            # Equivalence of default and keyword constructors
            @test PREMPolyModel(R, n, r, Vp', Vs', Rho', true, Vph', Vpv',
                Vsh', Vsv', Eta', true, qμ', qκ', fref) ==
                PREMPolyModel(r=r, vp=Vp, vs=Vs, density=Rho, vph=Vph,
                    vpv=Vpv, vsh=Vsh, vsv=Vsv, Qμ=qμ, Qκ=qκ, eta=Eta, fref=fref)
            @test PREMPolyModel(R, n, r, Vp', Vs', null, false, null,
                                     null, null, null, null, true, qμ', qκ', NaN) ==
                PREMPolyModel(r=r, vp=Vp, vs=Vs, Qμ=qμ, Qκ=qκ)
            @test PREMPolyModel(R, n, r, Vp', Vs', Rho', true, Vph', Vpv',
                                Vsh', Vsv', Eta', false, null, null, fref) ==
                PREMPolyModel(r=r, vp=Vp, vs=Vs, density=Rho, vph=Vph,
                              vpv=Vpv, vsh=Vsh, vsv=Vsv, eta=Eta, fref=fref)
            # ASCII arguments are the same and take precedence
            @test PREMPolyModel(r=r, Qmu=qμ, Qkappa=qκ) ==
                PREMPolyModel(r=r, Qμ=qμ, Qκ=qκ) ==
                PREMPolyModel(r=r, Qμ=Vs, Qκ=Vp, Qmu=qμ, Qkappa=qκ)
            # Arrays are converted to the correct shape
            @test_throws ArgumentError PREMPolyModel(r=r, vp=rand(3, n+1))
            Vp′ = rand(3, n)
            @test PREMPolyModel(r=r, vp=Vp′) == PREMPolyModel(R, n, r, Vp′,
                null, null, false, null, null, null, null, null, false, null, null, NaN)
        end
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
