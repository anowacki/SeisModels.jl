"""
# Moon

This module defines seismic models for the Earth's moon.
"""
module Moon

using ..SeisModels

"""
	ALL_MODELS

List of Earth model instances as `Symbol`s.

If adding new inbuilt Earth models, append the name of the new model to this list
"""
const ALL_MODELS = (:MOON_WEBER_2011,)

include("weber_2011.jl")

end # module
