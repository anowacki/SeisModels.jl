"""
# Earth

This module defines seismic models for the Earth.
"""
module Earth

using ..EarthModels

export AK135, PREM

include("ak135.jl")
include("prem.jl")

end # module
