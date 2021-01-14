#================
1D seismic models
================#
"""
Abstract type for radially-symmetric models of a planet.  1D seismic models should
be subtypes of this.
"""
abstract type SeisModel1D <: SeisModel end

const model_variables_SeisModel1D = (:vp, :vs, :density, :vph, :vpv, :vsh, :vsv,
    :eta, :Qμ, :Qκ)
const model_names_SeisModel1D = ("isotropic average Vp", "isotropic average Vs",
    "density", "Vph", "Vpv", "Vsh", "Vsv", "η", "Qμ", "Qκ")
const model_units_SeisModel1D = (" (km/s)", " (km/s)", " (g/cm^3)", " (km/s)",
    " (km/s)", " (km/s)", " (km/s)", "", "", "")

"""
    findlayer(m::SeisModel1D, r) -> layer_index

Return the layer index in which radius `r` km is located in the SeisModel1D `m`.
"""
function findlayer(m::SeisModel1D, r::Real)
    r < 0. && throw(DomainError(r, "radius cannot be negative"))
    r > m.a && throw(DomainError(r, "radius is greater than Earth radius for model ($(m.a) km)"))
    for i in 1:length(m.r)
        @inbounds r < m.r[i] && return i
    end
    length(m.r)
end

# Implementations of vp, vs, etc., evaluate here
include("prempolymodel.jl")
include("steppedlayeredmodel.jl")
include("linearlayeredmodel.jl")


## Derived quantities
"""
    bulk_modulus(m::SeisModel1D, r; depth=false) -> K

Return the bulk modulus `K` in GPa at radius `r` km in the model `m`.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function bulk_modulus(m::SeisModel1D, r; depth::Bool=false)
    depth && (r = radius(m, r))
    density(m, r)*vp(m, r)^2 - 4/3*shear_modulus(m, r)
end

"Newton's gravitational constant G = 6.67428e-11"
const NewtonG = 6.67428e-11

"""
    gravity(m::SeisModel1D, r; depth=false) -> g

Return the acceleration due to gravity, `g`, in m/s^2 at radius `r` in km.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function gravity(m::SeisModel1D, r; depth=false)
    depth && (r = radius(m, r))
    (r == 0) ? 0. : NewtonG*mass(m, r)/(r*1.e3)^2
end

# Individual implementations of mass are in model types' respective files
"""
    mass(m::SeisModel1D, r; depth=false) -> mass

Return the mass in kg between the centre of the model and the radius `r` km.

If `depth` is `true`, `r` is treated as a depth in km instead.

    mass(m) -> mass

Return the mass for the whole body in kg (between the centre of the model
and the surface).
"""
mass(m::SeisModel1D) = mass(m, surface_radius(m))

"""
    pressure(m::SeisModel1D, r; depth=false) -> p

Return the pressure `p` in Pa at radius `r` km.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function pressure(m::SeisModel1D, r; depth::Bool=false)
    depth && (r = radius(m, r))
    f(r) = pressure_integration_func(m, r)
    1.e3*quadgk(f, r, m.a)[1]
end
pressure_integration_func(m::SeisModel1D, r) = 1.e3*density(m, r)*gravity(m, r)

"""
    poissons_ratio(m, r; depth=false) -> ν

Return the Poisson's ratio for the model `m` given a radius `r` in km.
The calculation uses the isotropic average velocities at `r`.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function poissons_ratio(m::SeisModel1D, r; depth::Bool=false)
    depth && (r = radius(m, r))
    G = shear_modulus(m, r)
    K = bulk_modulus(m, r)
    (3K - 2G)/(6K + 2G)
end

"""
    shear_modulus(m::SeisModel1D, r; depth=false) -> μ

Return the shear modulus `μ` (often also called G) in GPa at radius `r` in the
model `m`.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function shear_modulus(m::SeisModel1D, r; depth::Bool=false)
    depth && (r = radius(m, r))
    vs(m, r)^2*density(m, r)
end

"""
    surface_mass(m::SeisModel1D, r; depth=false) -> mass

Return the mass in kg betwen radius `r` km and the surface.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function surface_mass(m::SeisModel1D, r; depth::Bool=false)
    depth && (r = radius(m, r))
    mass(m, m.a) - mass(m, r)
end

"""
    youngs_modulus(m, r; depth=false) -> E

Return the Young's modulus `E` in GPa for the model `m` given a radius `r` in km.

If `depth` is `true`, `r` is treated as a depth in km instead.
"""
function youngs_modulus(m::SeisModel1D, r; depth::Bool=false)
    depth && (r = radius(m, r))
    2*shear_modulus(m, r)*(1 + poissons_ratio(m, r))
end

"""
    moment_of_inertia(m, r0=0, r1=surface_radius(m)) -> I

Return the moment of interia `I` in kg m² for the model `m` between radii `r0`
and `r1` in km.
"""
moment_of_inertia(m::SeisModel1D, r0=0, r1=surface_radius(m)) =
    8/3*π*quadgk(r->1e3*density(m, r/1e3)*r^4, r0*1e3, r1*1e3)[1]
