# EarthModels.jl

## What is EarthModels.jl?
A [Julia](http://julialang.org) package for dealing with models of the Earth's
seismic properties.

Currently, only three kinds of one-dimensional models are supported.

Built in models are PREM and AK135.


## How to install
Although not registered as an official package, EarthModels.jl can be added to your
Julia install like so:

```julia
Pkg.clone("https://github.com/anowacki/EartModels.jl")
```

You then need only do

```julia
using EarthModels
```

and if that works, you're ready to go.


## How to use
### Model types
The module defines the EarthModel type and subtypes of this specify the kind of
model (i.e., symmetry, nature of basis function, etc.).

```julia
julia> using EarthModels

julia> subtypes(EarthModel)
1-element Array{Any,1}:
 EarthModels.EarthModel1D

julia> subtypes(EarthModel1D)
3-element Array{Any,1}:
 EarthModels.LinearLayeredModel 
 EarthModels.PREMPolyModel      
 EarthModels.SteppedLayeredModel
```

So, there are currently three types of models implemented, all 1D models, with
polynomial, linear or constant basis within each layer.

### Calculating properties

You can either create your own models by creating a new instance of one of the
immutable types, or use the inbuilt models.  For instance, for PREM, once can
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

julia> rho(AK135, AK135.a-20)
2.449
```

In the last example, we used the `a` field of the AK135 model type, which is the
radius of the Earth in km in the model, to calculate the density at 20 km depth.


## Getting help
Types and methods are documented, so at the REPL type `?` to get a `help?>`
prompt, and type the name of the function:

```julia
help?> PREMPolyModel
search: PREMPolyModel

  PREMPolyModel <: EarthModel1D

  Type describing the Earth as a set of layers within which properties vary according to a set of
  polynomials.

  Physical parameters are represented by arrays of size (order+1,n), where n is the number of
  layers, and order is the order of polynomial which is used to represent the parameter. Hence a
  constant layer has size (1,n), and to compute the value of x in layer i, for an Earth radius of
  a km, at a radius of r km, the expression is:

  val_x = x[i,1] + (r/a)*x[i,2] + (r/a)^2*x[i,3] ... (r/a)^order*x[i,order+1]

```

### Documentation
Documentation is a work in progress, but all useful commands are documented.
To see the list of commands, check the code, or in the REPL type `SAC.` then
press tab a couple of times to see all the module methods and variables.
Calling up the interactive help with give a useful description of each.


## Dependencies
None.  The package has been tested on Julia v0.4 and v0.5
