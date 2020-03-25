# PREMPolyModels

"""
    PREMPolyModel <: SeisModel1D

Type describing the Earth as a set of layers within which properties vary
according to a set of polynomials.

Physical parameters are represented by arrays of size `(order+1, n)`, where `n`
is the number of layers, and `order` is the order of polynomial which is
used to represent the parameter.  Hence a constant layer has size `(1,n)`,
and to compute the value of `x` in layer `i`, for an Earth radius of `a` km,
at a radius of `r` km, the expression is:

    val_x = x[i,1] + (r/a)*x[i,2] + (r/a)^2*x[i,3] ... (r/a)^order*x[i,order+1]

The reference frequency in Hz for the velocities in case of attenuation is
given by the `fref` field.  If `fref` is `NaN`, then no reference frequency
is defined for the model.
"""
struct PREMPolyModel <: SeisModel1D
    "Earth radius in km"
    a :: Float64
    "Number of distinct layers"
    n :: Int
    "Radii in km of top of each layer"
    r :: Vector{Float64}
    "Physical parameters"
    vp :: Array{Float64,2}
    vs :: Array{Float64,2}
    density :: Array{Float64,2}
    "Optional anisotropic parameters"
    aniso :: Bool
    vph :: Array{Float64,2}
    vpv :: Array{Float64,2}
    vsh :: Array{Float64,2}
    vsv :: Array{Float64,2}
    eta :: Array{Float64,2}
    "Optional attenuation parameters"
    attenuation :: Bool
    Qμ :: Array{Float64,2}
    Qκ :: Array{Float64,2}
    "Reference frequency in Hz"
    fref :: Float64
end

"""
    PREMPolyModel(; kwargs...) -> m

Construct a `PREMPolyModel` by specifying keyword arguments.  The length of
arrays defines the number of layers and the presence of arguments determines
whether the model has density, anisotropy and attenuation.

Note that radii are given in km, velocities in km/s, densities in g/cm³
and frequency in Hz.

### Keyword arguments
#### Radii
- `r`: The upper radii of each layer in km.  This is the only mandatory argument.
#### Isotropic properties
These must be given each as an `(norder + 1)`×`nlayers` `Matrix` of polynomial
coefficients.  Each parameter (`vp`, `density`, `Qμ`, etc.) may have a
different polynomial order, but all must have the correct number of layers
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
- `fref`: Reference Q frequency
"""
function PREMPolyModel(; r, vp=nothing, vs=nothing, density=nothing,
                       vph=nothing, vpv=nothing, vsh=nothing, vsv=nothing,
                       eta=nothing, Qμ=nothing, Qmu=nothing,
                       Qκ=nothing, Qkappa=nothing, fref=nothing)
    n = length(r)
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
    vp = _set_coeffs_size(PREMPolyModel, n, vp)
    vs = _set_coeffs_size(PREMPolyModel, n, vs)
    density = _set_coeffs_size(PREMPolyModel, n, density)
    vph = _set_coeffs_size(PREMPolyModel, n, vph)
    vsh = _set_coeffs_size(PREMPolyModel, n, vsh)
    vpv = _set_coeffs_size(PREMPolyModel, n, vpv)
    vsv = _set_coeffs_size(PREMPolyModel, n, vsv)
    eta = _set_coeffs_size(PREMPolyModel, n, eta)
    Qμ = _set_coeffs_size(PREMPolyModel, n, Qμ)
    Qκ = _set_coeffs_size(PREMPolyModel, n, Qκ)
    fref === nothing && (fref = NaN)
    PREMPolyModel(R, n, r, vp, vs, density, aniso, vph, vpv, vsh, vsv, eta,
        attenuation, Qμ, Qκ, fref)
end

# Check dimensions of input and resize if needed
_set_coeffs_size(::Type{PREMPolyModel}, n, v::Nothing) = Array{Float64}(undef, 0, 0)
function _set_coeffs_size(::Type{PREMPolyModel}, n, v::AbstractVector)
    length(v) == n && return permutedims(v)
    throw(ArgumentError("vector of coefficients must have length $n"))
end
function _set_coeffs_size(::Type{PREMPolyModel}, n, v::AbstractArray{A,2} where A)
    size(v, 2) == n ||
        throw(ArgumentError("second dimension of coefficients array must be $n"))
    v
end
_set_coeffs_size(::Type{PREMPolyModel}, n, v::A) where A =
    throw(ArgumentError("unexpected type of coefficients: $A"))


"""
    _evalpoly(x, coeffs, i)

Evaluate a polynomial whose coefficients are given in `coeffs[:,i]`
at `x`.  Coefficients should be orderd `a₀, a₁, ..., aₙ`.
"""
function _evalpoly(x, coeffs::AbstractArray{<:Real,2}, i)
    val = coeffs[1,i]
    for k in 2:size(coeffs, 1)
        @inbounds val = val + coeffs[k,i]*x^(k-1)
    end
    val
end

# Evaluation routines--all documented here
# TODO: Use Horner's method à la Base.@evalpoly
for (sym, name, unit) in zip(model_variables_SeisModel1D, model_names_SeisModel1D, model_units_SeisModel1D)
    # Velocities, which can be corrected for period if we have attenuation
    # and a reference frequency
    if sym in (:vp, :vs, :vph, :vpv, :vsh, :vsv)
        VP, VS = if sym in (:vp, :vs)
            :vp, :vs
        elseif sym in (:vph, :vsh)
            :vph, :vsh
        else
            :vpv, :vsv
        end
        _correct = Symbol(:_correct_attenuation_, String(sym)[2])
        @eval begin
            """
                $(split(string($sym), ".")[end])(m::SeisModel1D, r; depth=false, freq=reffrequency(m)) -> $($name)

            Return the value of $($name)$($unit) for model `m` at radius `r` km.

            If `depth` is `true`, then `r` is given as a depth in km instead.

            If `freq` is given, then the velocity is corrected for attenuation.
            This requires the model has attenuation and a reference frequency.
            """
            function ($sym)(m::PREMPolyModel, r::Real; depth::Bool=false, freq=nothing)
                length(m.$sym) > 0 || throw(ArgumentError("$($name) not defined for model"))
                depth && (r = radius(m, r))
                ir = findlayer(m, r)
                x = r/m.a
                val = _evalpoly(x, m.$(sym), ir)
                freq === nothing && return val
                hasattenuation(m) ||
                    throw(ArgumentError("cannot correct a nonattenuating model " *
                                        "for attenuation"))
                hasreffrequency(m) ||
                    throw(ArgumentError("no reference frequency defined for model"))
                freq == reffrequency(m) && return val
                _vp = $(sym == VP ? :val : :(_evalpoly(x, m.$(VP), ir)))
                _vs = $(sym == VS ? :val : :(_evalpoly(x, m.$(VS), ir)))
                E = 4/3*(_vs/_vp)^2
                qκ = 1/_evalpoly(x, m.Qκ, ir)
                qμ = 1/_evalpoly(x, m.Qμ, ir)
                $(_correct)(val, freq, reffrequency(m), E, qμ, qκ)
            end
        end
    else
        @eval begin
            """
                $(split(string($sym), ".")[end])(m::SeisModel1D, r; depth=false) -> $($name)
            
            Return the value of $($name)$($unit) for model `m` at radius `r` km.
            
            If `depth` is `true`, then `r` is given as a depth in km instead.
            """
            function ($sym)(m::PREMPolyModel, r::Real; depth::Bool=false)
                length(m.$sym) > 0 || throw(ArgumentError("$($name) not defined for model"))
                depth && (r = radius(m, r))
                ir = findlayer(m, r)
                x = r/m.a
                val = _evalpoly(x, m.$(sym), ir)
                val
            end
        end
    end
end
# Alternative names
const ρ = density
const Qmu = Qμ
const Qkappa = Qκ

"""
    evaluate(m::SeisModel1D, field::Symbol, r; depth=false) -> vals

Evaluate the model `m` at radius `r` km for the different property/ies in `field`,
returning a scalar for scalar input, and an array for array input.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function evaluate(m::PREMPolyModel, field::Symbol, r; depth::Bool=false)
    y = getfield(m, field)
    length(y) > 0 || throw(ArgumentError("'$field' not defined for model"))
    depth && (r = radius(m, r))
    ir = findlayer(m, r)
    x = r/m.a
    val = _evalpoly(x, y, ir)
    val
end

_correct_attenuation_s(vs, freq, reffreq, E, qμ, qκ) = vs*(1 - log(reffreq/freq)/π*qμ)
_correct_attenuation_p(vp, freq, reffreq, E, qμ, qκ) = vp*(1 - log(reffreq/freq)/π*((1 - E)*qκ + E*qμ))

function mass(m::PREMPolyModel, r; depth::Bool=false)
    hasdensity(m) || throw(ArgumentError("model does not contain density"))
    depth && (r = radius(m, r))
    l = findlayer(m, r)
    r *= 1.e3 # SI
    M = 0.
    for i = 1:l
        if i == 1
            R0 = 0.
        else
            R0 = m.r[i-1]*1.e3
        end
        if i < l
            R = m.r[i]*1.e3
        else
            R = r
        end
        for k in 1:size(m.density, 1)
            rho = m.density[k,i]/m.a^(k-1)*1.e3^(2-k)/(k+2)
            M += rho*(R^(k+2) - R0^(k+2))
        end
    end
    4*π*M
end
