"""
# Earth

This module defines seismic models for the Earth.
"""
module Earth

using ..SeisModels

"""
	ALL_MODELS

List of Earth model instances as `Symbol`s.

If adding new inbuilt Earth models, append the name of the new model to this list
"""
const ALL_MODELS = (:AK135, :EK137, :IASP91, :PREM, :PREM_NOOCEAN, :STW105)

include("ak135.jl")
include("ek137.jl")
include("iasp91.jl")
include("prem.jl")
include("stw105.jl")

end # module
