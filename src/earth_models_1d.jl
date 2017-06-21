#==============
1D Earth models
==============#
@compat abstract type EarthModel1D <: EarthModel end

const model_variables_EarthModel1D = (:vp, :vs, :rho, :vph, :vpv, :vsh, :vsv, :eta, :Qμ, :Qκ)
const model_names_EarthModel1D = ("Vp", "Vs", "density", "Vph", "Vpv", "Vsh", "Vsv", "η", "Qμ", "Qκ")
const model_units_EarthModel1D = (" (km/s)", " (km/s)", " (g/cm^3)", " (km/s)", " (km/s)", " (km/s)", " (km/s)", "", "", "")

"""
    findlayer(m::EarthModel1D, r) -> layer_index

Return the layer index in which radius `r` km is located in the EarthModel1D `m`.
"""
function findlayer(m::EarthModel1D, r::Real)
    r < 0. && error("Radius cannot be negative")
    r > m.a && error("Radius $r km is greater than Earth radius for model ($(m.a) km)")
    for i in 1:length(m.r)
        r < m.r[i] && return i
    end
    length(m.r)
end


"""
    PREMPolyModel <: EarthModel1D

Type describing the Earth as a set of layers within which properties vary
according to a set of polynomials.

Physical parameters are represented by arrays of size (order+1,n), where `n`
is the number of layers, and `order` is the order of polynomial which is
used to represent the parameter.  Hence a constant layer has size (1,n),
and to compute the value of x in layer i, for an Earth radius of a km,
at a radius of r km, the expression is:

    val_x = x[i,1] + (r/a)*x[i,2] + (r/a)^2*x[i,3] ... (r/a)^order*x[i,order+1]
"""
immutable PREMPolyModel <: EarthModel1D
    "Earth radius in km"
    a :: Float64
    "Number of distinct layers"
    n :: Int
    "Radii in km of top of each layer"
    r :: Vector{Float64}
    "Physical parameters"
    vp :: Array{Float64}
    vs :: Array{Float64}
    rho :: Array{Float64}
    "Optional anisotropic parameters"
    aniso :: Bool
    vph :: Array{Float64}
    vpv :: Array{Float64}
    vsh :: Array{Float64}
    vsv :: Array{Float64}
    eta :: Array{Float64}
    "Optional attenuation parameters"
    attenuation :: Bool
    Qμ :: Array{Float64}
    Qκ :: Array{Float64}
end

# Evaluation routines--all documented here
for (sym, name, unit) in zip(model_variables_EarthModel1D, model_names_EarthModel1D, model_units_EarthModel1D)
    @eval begin
        @doc """
            $(split(string($sym), ".")[end])(m::EarthModel1D, r) -> $($name)
        
        Return the value of $($name)$($unit) for model `m` at radius `r` km.
        """ ->
        function ($sym)(m::PREMPolyModel, r::Real)
            length(m.$sym) > 0 || error("$($name) not defined for model")
            ir = findlayer(m, r)
            x = r/m.a
            val = m.$(sym)[1,ir]
            for k in 2:size(m.$sym, 1)
                val = val + m.$(sym)[k,ir]*x^(k-1)
            end
            val
        end
    end
end
# Alternative names
const ρ = rho
const Qmu = Qμ
const Qkappa = Qκ

# Allow evaluation of multiple points simultaneously
for sym in model_variables_EarthModel1D
    @eval $sym(m::EarthModel1D, r::AbstractArray) = [$sym(m, R) for R in r]
end



"""
    SteppedLayeredModel <: EarthModel1D

A `SteppedLayeredModel` contains `n` layers with maximum radius `r` km,
each with a constant velocity.
"""
immutable SteppedLayeredModel <: EarthModel1D
    "Earth radius in km"
    a :: Float64
    "Number of layers"
    n :: Int
    "Radii in km of top of each layer"
    r :: Vector{Float64}
    "Physical parameters"
    vp :: Vector{Float64}
    vs :: Vector{Float64}
    rho :: Vector{Float64}
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

for sym in model_variables_EarthModel1D
    @eval function ($sym)(m::SteppedLayeredModel, r::Real)
        length(m.$sym) > 0 || error("$(split(string($sym), ".")[end]) not defined for model")
        m.$(sym)[findlayer(m, r)]
    end
end


"""
    LinearLayeredModel <: EarthModel1D

A `LinearLayeredModel` contains `n` points at which velocities are defined, with
linear interpolation between them.  Hence there are `n - 1` layers.

Discontinuities are represented by two layers with the same radii.
"""
immutable LinearLayeredModel <: EarthModel1D
    "Earth radius in km"
    a :: Float64
    "Number of layers"
    n :: Int
    "Radii in km of top of each layer"
    r :: Vector{Float64}
    "Physical parameters"
    vp :: Vector{Float64}
    vs :: Vector{Float64}
    rho :: Vector{Float64}
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

function findlayer(m::LinearLayeredModel, r::Real)
    r < 0. && error("Radius cannot be negative")
    r > m.a && error("Radius $r km is greater than Earth radius for model ($(m.a) km)")
    for i in 1:(length(m.r) - 1)
        r < m.r[i+1] && return i
    end
    return length(m.r) - 1
end

for sym in model_variables_EarthModel1D
    @eval function ($sym)(m::LinearLayeredModel, r::Real)
        length(m.$sym) > 0 || error("$(split(string($sym), ".")[end]) not defined for model")
        ir = findlayer(m, r)
        r0 = m.r[ir]
        r1 = m.r[ir+1]
        (m.$sym[ir]*(r1 - r) + m.$sym[ir+1]*(r - r0))/(r1 - r0)
    end
end


## Derived quantities
"""
    mass(m::EarthModel1D, r) -> mass

Return the mass in kg between the centre of the model and the radius `r` km.
"""
function mass(m::PREMPolyModel, r)
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
        for k in 1:size(m.rho, 1)
            rho = m.rho[k,i]/m.a^(k-1)*1.e3^(2-k)/(k+2)
            M += rho*(R^(k+2) - R0^(k+2))
        end
	end
	4*π*M
end

function mass(m::SteppedLayeredModel, r)
    l = findlayer(m, r)
    r *= 1.e3
    M = 0.
    for i in 1:l
        R0 = i == 1 ? 0. : m.r[i-1]*1.e3
        R = i < l ? m.r[i]*1.e3 : r
        M += m.rho[i]*(R^3 - R0^3)
    end
    4/3*π*M
end

function mass(m::LinearLayeredModel, r)
    l = findlayer(m, r)
    r *= 1.e3
    M = 0.
    for i in 1:l
        m.r[i] == m.r[i+1] && continue
        R0 = m.r[i]*1.e3
        R = i == l ? r : m.r[i+1]*1.e3
        dρ_dr = 1.e3*(m.rho[i+1] - m.rho[i])/(m.r[i+1]*1.e3 - R0)
        ρ0 = 1.e3*m.rho[i] - dρ_dr*R0
        M += ρ0*(R^3 - R0^3)/3 + dρ_dr*(R^4 - R0^4)/4
    end
    4*π*M
end

"""
    surface_mass(m::EarthModel1D, r) -> mass

Return the mass in kg betwen radius `r` km and the surface.
"""
surface_mass(m::EarthModel1D, r) = mass(m.a) - mass(r)

const NewtonG = 6.67428e-11
"""
    gravity(m::EarthModel1D, r) -> g

Return the acceleration due to gravity, `g`, in m/s^2 at radius `r` km.
"""
gravity(m::EarthModel1D, r) = (r == 0) ? 0. : NewtonG*mass(m, r)/(r*1.e3)^2

"""
    pressure(m::EarthModel1D, r) -> p

Return the pressure `p` in Pa at radius `r` km.
"""
function pressure(m::EarthModel1D, r)
    f(r) = pressure_integration_func(m, r)
    1.e3*quadgk(f, r, m.a)[1]
end
pressure_integration_func(m::EarthModel1D, r) = 1.e3*rho(m, r)*gravity(m, r)

"""
    bulk_modulus(m::EarthModel1D, r) -> K

Return the bulk modulus `K` in GPa at radius r in the model `m`.
"""
bulk_modulus(m::EarthModel1D, r) = rho(m, r)*vp(m, r)^2 - 4/3*shear_modulus(m, r)
bulk_modulus(m::EarthModel1D, r::AbstractArray) = [bulk_modulus(m, R) for R in r]

"""
    shear_modulus(m::EarthModel1D, r) -> μ

Return the shear modulus `μ` in GPa at radius r in the model `m`.
"""
shear_modulus(m::EarthModel1D, r) = vs(m, r)^2*rho(m, r)
shear_modulus(m::EarthModel1D, r::AbstractArray) = [shear_modulus(m, R) for R in r]
