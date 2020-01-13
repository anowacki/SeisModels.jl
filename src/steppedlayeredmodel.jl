"""
    SteppedLayeredModel <: SeisModel1D

A `SteppedLayeredModel` contains `n` layers each with constant properties.
"""
struct SteppedLayeredModel <: SeisModel1D
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
    SteppedLayeredModel(; kwargs...) -> m

Construct a `SteppedLayeredModel` by specifying keyword arguments.  The length of
arrays defines the number of layers and the presence of arguments determines
whether the model has density, anisotropy and attenuation.

Note that radii are given in km, velocities in km/s, densities in g/cm³
and frequency in Hz.

### Keyword arguments
#### Radii
- `r`: The upper radius of a set of constant layers in km.  The uppermost
  layer radius is taken to be the planetary radius.  This is the only
  mandatory argument.
#### Isotropic properties
These must be given each as a `Vector` of values for each layer.
Each parameter (`vp`, `density`, `Qμ`, etc.) must have the correct number
of layers.
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
function SteppedLayeredModel(; r, vp=nothing, vs=nothing, density=nothing,
                             vph=nothing, vpv=nothing, vsh=nothing, vsv=nothing,
                             eta=nothing, Qμ=nothing, Qmu=nothing,
                             Qκ=nothing, Qkappa=nothing, fref=nothing)
    n = length(r)
    all(x->x>0, diff(r)) || throw(ArgumentError("radii must increase monotonically"))
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
    vp = _set_coeffs_size(SteppedLayeredModel, n, vp)
    vs = _set_coeffs_size(SteppedLayeredModel, n, vs)
    density = _set_coeffs_size(SteppedLayeredModel, n, density)
    vph = _set_coeffs_size(SteppedLayeredModel, n, vph)
    vsh = _set_coeffs_size(SteppedLayeredModel, n, vsh)
    vpv = _set_coeffs_size(SteppedLayeredModel, n, vpv)
    vsv = _set_coeffs_size(SteppedLayeredModel, n, vsv)
    eta = _set_coeffs_size(SteppedLayeredModel, n, eta)
    Qμ = _set_coeffs_size(SteppedLayeredModel, n, Qμ)
    Qκ = _set_coeffs_size(SteppedLayeredModel, n, Qκ)
    SteppedLayeredModel(R, n, r, vp, vs, density, aniso, vph, vpv,
        vsh, vsv, eta, attenuation, Qμ, Qκ)
end

_set_coeffs_size(::Type{SteppedLayeredModel}, n, v::Nothing) = []
function _set_coeffs_size(::Type{SteppedLayeredModel}, n, v::Vector)
    length(v) == n || throw(ArgumentError("number of value must be $n"))
    v
end

for sym in model_variables_SeisModel1D
    @eval function ($sym)(m::SteppedLayeredModel, r::Real; depth::Bool=false)
        length(m.$sym) > 0 ||
            throw(ArgumentError("'$split(string($sym), ".")[end]' not defined for model"))
        depth && (r = radius(m, r))
        m.$(sym)[findlayer(m, r)]
    end
end

function evaluate(m::SteppedLayeredModel, field::Symbol, r; depth::Bool=false)
    length(getfield(m, field)) > 0 ||
        throw(ArgumentError("$field' not defined for model"))
    depth && (r = radius(m, r))
    getfield(m, field)[findlayer(m, r)]
end

function mass(m::SteppedLayeredModel, r; depth::Bool=false)
    hasdensity(m) || throw(ArgumentError("model does not contain density"))
    depth && (r = radius(m, r))
    l = findlayer(m, r)
    r *= 1.e3
    M = 0.
    for i in 1:l
        R0 = i == 1 ? 0. : m.r[i-1]*1.e3
        R = i < l ? m.r[i]*1.e3 : r
        M += 1.e3*m.density[i]*(R^3 - R0^3)
    end
    4/3*π*M
end
