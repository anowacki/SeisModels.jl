module EarthModels

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
    g,
    mass,
    pressure,
    surface_mass,

    # Models
    AK135,
    PREM

# Type of which all others are subtypes
abstract EarthModel

## 1D models
include("earth_models_1d.jl")
# Predefined 1D models
include("prem.jl")
include("ak135.jl")

end # module
