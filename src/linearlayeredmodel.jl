"""
    LinearLayeredModel <: SeisModel1D

A `LinearLayeredModel` contains `n` points at which velocities are defined, with
linear interpolation between them.  Hence there are `n - 1` layers.

Discontinuities are represented by two layers with the same radii.
"""
struct LinearLayeredModel <: SeisModel1D
    "Earth radius in km"
    a :: Float64
    "Number of layers"
    n :: Int
    "Radii in km of top of each layer"
    r :: Vector{Float64}
    "Physical parameters"
    vp :: Vector{Float64}
    vs :: Vector{Float64}
    density :: Vector{Float64}
    "Optional anisotropic parameters"
    aniso :: Bool
    vph :: Vector{Float64}
    vpv :: Vector{Float64}
    vsh :: Vector{Float64}
    vsv :: Vector{Float64}
    eta :: Vector{Float64}
    "Optional attenuation parameters"
    attenuation :: Bool
    Qμ :: Vector{Float64}
    Qκ :: Vector{Float64}
end

"""
    LinearLayeredModel(; kwargs...) -> m

Construct a `LinearLayeredModel` by specifying keyword arguments.  The length of
arrays defines the number of layers and the presence of arguments determines
whether the model has density, anisotropy and attenuation.

Note that radii are given in km, velocities in km/s, densities in g/cm³
and frequency in Hz.

### Keyword arguments
#### Radii
- `r`: The radius of a set of nodes ranging from the centre of the planet to
  the surface, in km.  When the same value of `r` is repeated twice, a
  discontinuity is introduced at that radius.  The first value of `r` must
  be 0 km.  This is the only mandatory argument.
#### Isotropic properties
These must be given each as a `Vector` of values at each radial node.
Each parameter (`vp`, `density`, `Qμ`, etc.) must have the correct number of layers.
- `vp`: P-wave velocity
- `vs`: S-wave velocity
- `density`: Material density
#### Anisotropic parameters
If the model is anisotropic and the radial anisotropy parameters are all
provided, then the Voigt average velocities are constructed automatically
unless `vp` and/or `vs` are also given separately.
- `vpv`: Vertical P-wave velocity
- `vph`: Horizontal P-wave velocity
- `vsv`: Vertical S-wave velocity
- `vsh`: Horizontally-polarised S-wave velocity
- `eta`: The ratio C₁₃/(Vph - 2Vsv), where cᵢⱼ is the Voigt elasticity
  matrix, the 1-direction is horizontal and the 3-direction is radial
#### Attenuation parameters
ASCII versions (`Qmu`, `Qkappa`) of these exist for ease of use as well as
the Unicode versions (`Qμ`, `Qκ`).  If both are supplied, an error is thrown.
- `Qμ` or `Qmu`: Shear modulus quality factor
- `Qκ` or `Qkappa`: Bulk modulus quality factor
"""
function LinearLayeredModel(; r, vp=nothing, vs=nothing, density=nothing,
                            vph=nothing, vpv=nothing, vsh=nothing, vsv=nothing,
                            eta=nothing, Qμ=nothing, Qmu=nothing,
                            Qκ=nothing, Qkappa=nothing, fref=nothing)
    n = length(r)
    first(r) == 0 || throw(ArgumentError("first value of r must be 0 km"))
    all(x->x>=0, diff(r)) || throw(ArgumentError("radii must increase monotonically"))
    R = last(r)
    aniso = attenuation = false
    # Anisotropy
    if all(!isnothing, (vph, vpv, vsh, vsv, eta))
        aniso = true
        # Set isotropic velocities automatically if not provided
        if any(isnothing, (vp, vs))
            voigt_vels = voigt_velocities.(vpv, vsv, vph, vsh, eta)
            vp === nothing && (vp = first.(voigt_vels))
            vs === nothing && (vs = last.(voigt_vels))
        end
    end
    # Attenuation
    # Default to ASCII versions if given
    Qmu !== nothing && (Qμ = Qmu)
    Qkappa !== nothing && (Qκ = Qkappa)
    if all(!isnothing, (Qμ, Qκ))
        attenuation = true
    end
    # Convert to correct shape, and check for length
    vp = _set_coeffs_size(LinearLayeredModel, n, vp)
    vs = _set_coeffs_size(LinearLayeredModel, n, vs)
    density = _set_coeffs_size(LinearLayeredModel, n, density)
    vph = _set_coeffs_size(LinearLayeredModel, n, vph)
    vsh = _set_coeffs_size(LinearLayeredModel, n, vsh)
    vpv = _set_coeffs_size(LinearLayeredModel, n, vpv)
    vsv = _set_coeffs_size(LinearLayeredModel, n, vsv)
    eta = _set_coeffs_size(LinearLayeredModel, n, eta)
    Qμ = _set_coeffs_size(LinearLayeredModel, n, Qμ)
    Qκ = _set_coeffs_size(LinearLayeredModel, n, Qκ)
    LinearLayeredModel(R, n, r, vp, vs, density, aniso, vph, vpv,
        vsh, vsv, eta, attenuation, Qμ, Qκ)
end

_set_coeffs_size(::Type{LinearLayeredModel}, n, v::Nothing) = []
function _set_coeffs_size(::Type{LinearLayeredModel}, n, v::Vector)
    length(v) == n || throw(ArgumentError("number of values must be $n"))
    v
end

function findlayer(m::LinearLayeredModel, r::Real)
    r < 0. && error("Radius cannot be negative")
    r > m.a && error("Radius $r km is greater than Earth radius for model ($(m.a) km)")
    for i in 1:(length(m.r) - 1)
        r < m.r[i+1] && return i
    end
    return length(m.r) - 1
end

for sym in model_variables_SeisModel1D
    @eval function ($sym)(m::LinearLayeredModel, r::Real; depth::Bool=false)
        length(m.$sym) > 0 ||
            throw(ArgumentError("'$(split(string($sym), ".")[end])' not defined for model"))
        depth && (r = radius(m, r))
        ir = findlayer(m, r)
        r0 = m.r[ir]
        r1 = m.r[ir+1]
        (m.$sym[ir]*(r1 - r) + m.$sym[ir+1]*(r - r0))/(r1 - r0)
    end
end

function evaluate(m::LinearLayeredModel, field::Symbol, r; depth::Bool=false)
    y = getfield(m, field)
    length(y) > 0 ||
        throw(ArgumentError("'$field' not defined for model"))
    depth && (r = radius(m, r))
    ir = findlayer(m, r)
    r0 = m.r[ir]
    r1 = m.r[ir+1]
    (y[ir]*(r1 - r) + y[ir+1]*(r - r0))/(r1 - r0)
end

function mass(m::LinearLayeredModel, r; depth::Bool=false)
    hasdensity(m) || throw(ArgumentError("model does not contain density"))
    depth && (r = radius(m, r))
    l = findlayer(m, r)
    r *= 1.e3
    M = 0.
    for i in 1:l
        m.r[i] == m.r[i+1] && continue
        R0 = m.r[i]*1.e3
        R = i == l ? r : m.r[i+1]*1.e3
        dρ_dr = 1.e3*(m.density[i+1] - m.density[i])/(m.r[i+1]*1.e3 - R0)
        ρ0 = 1.e3*m.density[i] - dρ_dr*R0
        M += ρ0*(R^3 - R0^3)/3 + dρ_dr*(R^4 - R0^4)/4
    end
    4*π*M
end
