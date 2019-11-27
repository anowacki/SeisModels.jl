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
	hasdensity(m) -> ::Bool

Return `true` if density is defined for model `m`, and `false` otherwise.
"""
hasdensity(m::SeisModel) = !isempty(m.density)

"""
    hasreffrequency(m::PREMPolyModel) -> ::Bool

Return `true` if a reference frequency is defined for `PREMPolyModel` `m`,
and `false` otherwise.
"""
hasreffrequency(m::PREMPolyModel) = !isnan(m.fref)

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
    reffrequency(m::PREMPolyModel) = f

Return the reference frequency for a `PREMPolyModel` `m` in Hz.
"""
function reffrequency(m::PREMPolyModel)
    hasreffrequency(m) ||
        throw(ArgumentError("model does not have a reference frequency"))
    m.fref
end

"""
    surface_radius(m) -> radius

Return the surface radius in km of the `SeisModel` `m`.
"""
surface_radius(m::SeisModel) = m.a
