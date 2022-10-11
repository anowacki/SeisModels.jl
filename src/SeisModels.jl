"""
# SeisModels

Give various physical properties of the interior of the Earth for a range of
models.

### Units

Throughout, the module uses the following unit conventions:

- Distances are in km
- Velocities are in km/s
- Densities are in g/cm³
- Accelerations are in m/s²
- Pressures (from `pressure`) are in Pa
- Moduli (e.g., from `youngs_modulus`) are in GPa
- Masses in are in kg

### Physical properties

All models should have:

- P-wave velocity (accessed with the `vp` function) in km/s
- S-wave velocity (`vs`) in km/s
- (Average) Earth surface radius (`surface_radius`) in km

Models may also provide:

- Density (`density`) in g/cm^3
- Radial anisotropy:
  - VPV, VPH, VSV and VSH for radially-anistropic models (`vpv`, etc.) in km/s
  - Radial anisotropy parameter η (`eta`)
- Attenuation:
  - Shear wave quality factor (`Qμ`)
  - Bulk sound quality factor (`Qκ`)

(For anisotropic models, the values of VP and VS rather than VPV etc. should
be some average—perhaps the Voigt average—of the anisotropic velocities.)

Derived properties for all models:

- Radius from centre of Earth to a given depth (`depth`) in km
- Bulk and shear moduli (`bulk_modulus`, `shear_modulus`) in GPa
- Young's modulus (`youngs_modulus`) in GPa
- Poisson's ratio

Derived properties can be computed for all models with density:

- Mass up to radius r (`mass`) in kg
- Mass from radius r to the surface (`surface_mass`) in kg
- Pressure (`pressure`) in Pa
- Acceleration due to gravity g (`gravity`) in m/s^2

### Inbuilt models

- Earth
  - `AK135`
  - `EK137`
  - `IASP91`
  - `PREM`
  - `PREM_NOOCEAN`
  - `STW105`
- Moon
  - `MOON_WEBER_2011`
"""
module SeisModels

using DelimitedFiles: readdlm, writedlm
using Printf: @printf

using QuadGK: quadgk

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
    hasreffrequency,
    isanisotropic,
    radius,
    reffrequency,
    surface_radius,
    discontinuities,

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
    write_mineos,
    read_tvel,
    write_tvel

@static if VERSION < v"1.1"
    isnothing(x) = x === nothing
end

"""
Abstract supertype of all models of the Earth in the `SeisModels` module.
All models should be a subtype of this or one of `SeisModel`'s abstract
subtypes.
"""
abstract type SeisModel end

# Allow models to be broadcasted as scalars
Base.broadcastable(x::SeisModel) = Ref(x)

# Default equality and approximate equality comparisons
Base.:(==)(m1::T, m2::T) where {T<:SeisModel} =
    all(isequal(getfield(m1, f), getfield(m2, f)) for f in fieldnames(T))
Base.:(≈)(m1::T, m2::T; kwargs...) where {T<:SeisModel} =
    all(≈(getfield(m1, f), getfield(m2, f); kwargs...) for f in fieldnames(T))

## 1D models
include("models_1d.jl")

## Basic properties of models
include("basic_properties.jl")

## Conversion between models
include("conversion.jl")

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

end # module
