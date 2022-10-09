# SeisModels.jl

## Build status
[![Build Status](https://github.com/anowacki/SeisModels.jl/workflows/CI/badge.svg)](https://github.com/anowacki/SeisModels.jl/actions)
[![Code coverage](https://codecov.io/gh/anowacki/SeisModels.jl/branch/master/graph/badge.svg?token=OX6WPKMG7G)](https://codecov.io/gh/anowacki/SeisModels.jl)

## Documentation
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://anowacki.github.io/SeisModels.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://anowacki.github.io/SeisModels.jl/dev)

## Publication status
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02043/status.svg)](https://doi.org/10.21105/joss.02043)


## What is SeisModels.jl?
A [Julia](http://julialang.org) package for dealing with models of the Earth's
(and other quasi-1D planets') seismic properties.

Currently, only three kinds of one-dimensional models are supported, but all model
parameterisations and models are acceptable for inclusion.  Contributions
are welcome.

Built in models are [listed in the documentation](https://anowacki.github.io/SeisModels.jl/dev/inbuilt_models/).


## How to install
SeisModels.jl can be added to your Julia install like so:

```julia
julia> import Pkg; Pkg.add("SeisModels")
```


## How to use
### Model types
The module defines the SeisModel type and subtypes of this specify the kind of
model (i.e., symmetry, nature of basis function, etc.).

```julia
julia> using SeisModels

julia> subtypes(SeisModel)
1-element Array{Any,1}:
 SeisModel1D

julia> subtypes(SeisModel1D)
3-element Array{Any,1}:
 LinearLayeredModel
 PREMPolyModel
 SteppedLayeredModel
```

So, there are currently three types of models implemented, all 1D models, with
polynomial, linear or constant basis within each layer.

### Calculating properties

You can either create your own models by creating a new instance of one of the
immutable types, or use the inbuilt models.  For instance, for PREM, one can
evaluate at an arbitrary radius:

* *V*<sub>P</sub>
* *V*<sub>S</sub>
* density *&rho;*
* anisotropic parameters *V*<sub>PH</sub>, *V*<sub>PV</sub>, *V*<sub>SH</sub>,
  *V*<sub>SV</sub> and *&eta;*
  
Calculate these by calling the function with the model as the first argument:

```julia
julia> vp(PREM, 3500)
13.71171655163979

julia> Qκ(PREM, 1000)
1327.7

julia> density(AK135, radius(AK135, 20))
2.449
```

In the last example, we used the `radius` function to convert depth in the AK135 model
to radius and calculate the density at 20 km depth.  Some functions also accept the
`depth` keyword argument to instead evaluate properties at a point below the surface:

```julia
julia> density(AK135, radius(AK135, 20)) == density(AK135, 20, depth=true)
true
```

You can also evaluate values programmatically (i.e., where the parameter of
interest is a variable) by using the exported `evaluate` function, and broadcast
the call to get multiple values:

```julia
julia> evaluate(AK135, :vp, 3580)
13.653094354838709

julia> parameters = (:vp, :vs, :density);

julia> evaluate.(AK135, parameters, 3680)
(13.591187999999999, 7.226264, 5.4003499999999995)
```

### Model input and output
Support for reading and writing model files is currently limited.  However, SeisModels
does support reading and writing of
[Mineos](https://geodynamics.org/cig/software/mineos/)-format &lsquo;tabular&rsquo; models
(i.e., `SteppedLayeredModel`s) via the `read_mineos` and `write_mineos` functions.


## Reference
### Exported types
- `SeisModel`: Abstract supertype of all models
  - `SeisModel1D`: Abstract supertype of 1D models
    - `LinearLayeredModel`: 1D model with linearly-varying properties between node points
    - `PREMPolyModel`: 1D model defined by PREM-style polynomials (of arbitrary degree)
    - `SteppedLayeredModel`: 1D model with constant properties between node points

### Exported model instances
- Earth
  - `AK135`
  - `IASP91`
  - `PREM`
  - `PREM_NOOCEAN`
  - `STW105`
- Moon
  - `MOON_WEBER_2011`

### Exported functions
#### Model properties
- `depth`: Return depth in km given a radius and model
- `hasattenuation`: Whether a model includes attenuation
- `hasdensity`: Whether a model includes density
- `isanisotropic`: Whether a model is anisotropic
- `radius`: Return radius in km given a depth and model
- `surface_radius`: Radius in km of planet

#### Evaluation functions
- `evaluate`: Evaluate a given field for a model at any radius
- `vp`: P-wave velocity in km/s
- `vs`: S-wave velocity in km/s
- `density`: Density in g/cm^3
- `vph`: Horizontal P-wave velocity in km/s
- `vpv`: Vertical (radial) P-wave velocity in km/s
- `vsh`: Horizontally-polarised S-wave velocity in km/s
- `vsv`: Vertically-polarised S-wave velocity in km/s
- `eta`: Anisotropic parameter
- `Qμ`, `Qmu`: Shear quality factor
- `Qκ`, `Qkappa`: Bulk quality factor

#### Derived properties
- `bulk_modulus`: Bulk modulus (_K_) in Pa
- `gravity`: Acceleration due to gravity in m/s^2 at a given radius
- `mass`: Mass in kg from centre of model to a given radius
- `moment_of_inertia`: MOI in kg m^2
- `poissons_ratio`: Poisson's ratio
- `pressure`: Pressure in Pa
- `shear_modulus`: Shear modulus (_G_) in Pa
- `surface_mass`: Mass between two radii
- `youngs_modulus`: Young's modulus in Pa

#### IO
- `read_mineos`: Read Mineos tabular-format file
- `write_mineos`: Write Mineos tabular-format file
- `read_tvel`: Write tvel-format file
- `write_tvel`: Write tvel-format file


## Getting help
Types and methods are documented, so at the REPL type `?` to get a `help?>`
prompt, and type the name of the function:

```julia
help?> PREMPolyModel
search: PREMPolyModel

  PREMPolyModel <: SeisModel1D

  Type describing the Earth as a set of layers within which properties vary according to a set of
  polynomials.

  Physical parameters are represented by arrays of size (order+1,n), where n is the number of
  layers, and order is the order of polynomial which is used to represent the parameter. Hence a
  constant layer has size (1,n), and to compute the value of x in layer i, for an Earth radius of
  a km, at a radius of r km, the expression is:

  val_x = x[i,1] + (r/a)*x[i,2] + (r/a)^2*x[i,3] ... (r/a)^order*x[i,order+1]

```

## Citing
If you use SeisModels.jl for your work, please cite the following paper:
- Nowacki, A., 2020. SeisModels.jl: A Julia package for models of the
  Earth’s interior. Journal of Open Source Software 5, 2043.
  doi:[10.21105/joss.02043](https://doi.org/10.21105/joss.02043)
