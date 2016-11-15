module EarthModels

import Base: eta

export
    # Model types
    EarthModel,
    EarthModel1D,
    LinearLayeredModel,
    PREMPolyModel,
    SteppedLayeredModel,
    
    # Evluation functions
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
    
    # Models
    AK135,
    PREM


abstract EarthModel

include("earth_models_1d.jl")

#================
Predefined models
================#
include("prem.jl")
include("ak135.jl")

end # module
