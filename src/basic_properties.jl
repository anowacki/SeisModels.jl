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

"""
    voigt_kappa_mu(A, C, L, N, F) -> κ, μ

Return the Voigt average bulk modulus `κ` and shear modulus `μ` given the Love parameters
of a transversely isotropic medium.  If A, C, etc., contain density, then κ and μ will be
true moduli, whilst if the Love parameters are provided as velocities, then these
are density-normalised 'moduli' with dimensions of L²⋅T².
"""
voigt_kappa_mu(A, C, L, N, F) = (4A + C + 4F - 4N)/9, (A + C - 2F + 5N + 6L)/15

"""
    voigt_velocities(vpv, vsv, vph, vsh, η) -> vp, vs

Return the Voigt average `vp` and `vs`, in the same units as the input, given four
velocities and the radial anisotropy parameter `η`.
"""
function voigt_velocities(vpv, vsv, vph, vsh, η)
    A = vph^2
    C = vpv^2
    L = vsv^2
    N = vsh^2
    F = η*(A - 2L)
    κ, μ = voigt_kappa_mu(A, C, L, N, F) # In km²/s²
    vp = sqrt(κ + 4/3*μ)
    vs = sqrt(μ)
    vp, vs
end

"""
    discontinuities(m::LinearLayeredModel; depths=false) -> radii, indices

Return the `radii` (or depths if `depths=true`; both in km) and `indices`
of velocity discontinuities for the 1D `LinearLayeredModel` `m`.

Discontinuities are marked in `LinearLayeredModel`s by two repeated radii,
and the index returned for each one is the lower index.  Therefore
`indices` gives the properties of the model below the discontinuities,
and `indices .+ 1` gives the properties above.

# Example
```
julia> radii, inds = discontinuities(AK135)
([1217.5, 3479.5, 3631.0, 5711.0, 5961.0, 6161.0, 6336.0, 6351.0], [24, 69, 73, 116, 122, 127, 132, 134])

julia> index_icb = inds[1]
24

julia> AK135.vp[index_icb], AK135.vp[index_icb + 1]
(11.0427, 10.289)
```
"""
function discontinuities(m::LinearLayeredModel; depths=false)
    # Discontinuities are represented as duplicated radii in LinearLayeredModels.
    duplicates = Float64[]
    indices = Int[]
    for i = 2:length(m.r)
        if (
            isequal(m.r[i], m.r[i-1]) &&
            (length(duplicates) == 0 || !isequal(duplicates[end], m.r[i]))
        )
            push!(duplicates, m.r[i])
            push!(indices, i-1)
        end
    end
    depths ? (reverse(depth.(m, duplicates)), reverse(indices)) : (duplicates, indices)
end
