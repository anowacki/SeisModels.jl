"""
# SeisModels

Give various physical properties of the interior of the Earth for a range of
models.

### Physical properties

All models should have:

- P-wave velocity (accessed with the `vp` function) in km/s
- S-wave velocity (`vs`) in km/s
- Density (`density`) in g/cm^3
- (Average) Earth surface radius (`surface_radius`) in km

Models may also provide:

- VPV, VPH, VSV and VSH for radially-anistropic models (`vpv`, etc.) in km/s
- Radial anisotropy parameter η (`eta`)
- Shear wave quality factor (`Qμ`)
- Bulk sound quality factor (`Qκ`)

Derived properties can be computed for all models:

- Mass up to radius r (`mass`) in kg
- Mass from radius r to the surface (`surface_mass`) in kg
- Pressure (`pressure`) in Pa
- Acceleration due to gravity g (`gravity`) in m/s^2
- Bulk and shear module (`bulk_modulus`, `shear_modulus`) in GPa
- Radius from centre of Earth to a given depth (`depth`) in km

### Currently implemented models

- PREM
- IASP91
- AK135
"""
module SeisModels

import QuadGK: quadgk

export
    # Model types
    SeisModel,
    SeisModel1D,
    LinearLayeredModel,
    PREMPolyModel,
    SteppedLayeredModel,

    # Model properties
    depth,
    hasattenuation,
    hasdensity,
    isanisotropic,
    radius,
    surface_radius,

    # Evaluation functions
    evaluate,
    vp,
    vs,
    density,
    vph,
    vpv,
    vsh,
    vsv,
    eta,
    Qμ, Qmu,
    Qκ, Qkappa,

    # Derived properties
    bulk_modulus,
    gravity,
    mass,
    moment_of_inertia,
    poissons_ratio,
    pressure,
    shear_modulus,
    surface_mass,
    youngs_modulus,

    # IO
    read_mineos,
    write_mineos

"""
Abstract supertype of all models of the Earth in the `SeisModels` module.
All models should be a subtype of this or one of `SeisModel`'s abstract
subtypes.
"""
abstract type SeisModel end

# Allow models to be broadcasted as scalars
Base.broadcastable(x::SeisModel) = Ref(x)

## Basic properties of models
include("basic_properties.jl")

## 1D models
include("earth_models_1d.jl")

## IO
include("io.jl")

## Predefined models
# For each `planet`, models are defined in a module in "$(planet)/$(planet).jl".
# Each module must include a list of model symbols in ALL_MODELS.
for planet in (:Earth, :Moon)
    include(joinpath(@__DIR__, string(planet), string(planet, ".jl")))
    for model in getfield(SeisModels, planet).ALL_MODELS
        @eval begin
            using .$(planet): $(model)
            export $(model)
        end
    end
end

## Deprecated functions
rho(m::SeisModel, args...; kwargs...) =
    (@warn("`rho` is deprecated; use `density` instead"); density(m, args...; kwargs...))
export rho

end # module
