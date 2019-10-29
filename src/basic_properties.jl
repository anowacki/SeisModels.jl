# This file provides getter-type methods which essentially define the fields
# needed in any concrete subtype of SeisModel.

"""
    depth(m, radius) -> depth

Return the `depth` in km for the model `m` given a `radius` in km.
"""
depth(m::SeisModel, radius) = surface_radius(m) - radius

"""
    hasattenuation(m) -> ::Bool

Return `true` if the model `m` includes attenuation and `false` otherwise.
"""
hasattenuation(m::SeisModel) = m.attenuation

"""
    isanisotropic(m) -> ::Bool

Return `true` if the model `m` is anisotropic and `false` otherwise.
"""
isanisotropic(m::SeisModel) = m.aniso

"""
    radius(m, depth) -> radius

Return the `radius` in km for the model `m` given a `depth` in km.
"""
radius(m::SeisModel, depth) = surface_radius(m) - depth

"""
    surface_radius(m) -> radius

Return the surface radius in km of the `SeisModel` `m`.
"""
surface_radius(m::SeisModel) = m.a
