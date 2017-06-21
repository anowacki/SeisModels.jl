__precompile__()

module EarthModels

import Compat: @compat

import Base: eta

export
    # Model types
    EarthModel,
    EarthModel1D,
    LinearLayeredModel,
    PREMPolyModel,
    SteppedLayeredModel,

    # Evaluation functions
    vp,
    vs,
    rho,
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
    pressure,
    shear_modulus,
    surface_mass,

    # Models
    AK135,
    PREM

# Type of which all others are subtypes
@compat abstract type EarthModel end

## 1D models
include("earth_models_1d.jl")
# Predefined 1D models
include("prem.jl")
include("ak135.jl")

end # module
