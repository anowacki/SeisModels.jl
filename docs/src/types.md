# Model types

## Type hierarchy
The module defines the SeisModel type and subtypes of this specify the kind of
model (i.e., symmetry, nature of basis function, etc.).

The current type hierarchy is as follows:

- `SeisModel`: Abstract type of all seismic models
  - `SeisModel1D`: Abstract type of all one-dimensonal (radially-symmetric)
    seismic models
    - `SteppedLayeredModel`: Type defined by layers with contstant properties.
    - `LinearLayeredModel`: Type defined by layers with constant gradient.
    - `PREMPolyModel`: Type defined by polynomials.

As you can see, there are currently three types of models implemented, all
1D models, with polynomial, linear or constant basis within each layer.

### `SteppedLayeredModel`

A `SteppedLayeredModel` is composed of a series of layers with constant
properties in each layer.  The radii define the top of each layer.

An example of this is Weber et al.’s (2011) moon model:
```@eval
using SeisModels
using Plots: plot, plot!, vline!, savefig
p = plot(legend=:bottom, framestyle=:box, grid=false, title="SteppedLayeredModel example: the moon", xlabel="Radius (km)", ylabel="Velocity (km/s) or density (g/cm^3)")
for (f,l) in zip((vp, vs, SeisModels.density), ("Vp", "Vs", "Density"))
    plot!(p, r->f(MOON_WEBER_2011, r), 0, surface_radius(MOON_WEBER_2011), label=l, line=2)
end
vline!(p, MOON_WEBER_2011.r, label="Layer tops", l=(0.5, :dash))
savefig(p, "SteppedLayeredModel_example.svg")
nothing
```
![Example of a SteppedLayeredModel](SteppedLayeredModel_example.svg)

```@docs
SteppedLayeredModel
```

See the [constructor](@ref SeisModels.SteppedLayeredModel()) for details on creating
`SteppedLayeredModel`s.

### `LinearLayeredModel`

A `LinearLayeredModel` is described by radial nodes, between which
the properties vary linearly.  Discontinuities are added by setting
two subsequent nodes to have the same radii but different properties.

For example, a two-layer model for a planetesimal with radius 1 km
and a liquid core at 0.3 km, with linearly varying properties, could
be given by:

```@repl example
using SeisModels
m = LinearLayeredModel(r=[0, 0.3, 0.3, 1], vp=[1, 1, 1.8, 1], vs=[0, 0, 1, 0.7])
```

You can see that the nodes are used to define a constant lower layer
and linearly-varying upper layer:
```@eval
using SeisModels
using Plots: plot, plot!, scatter!, savefig
m = LinearLayeredModel(r=[0, 0.3, 0.3, 1], vp=[1, 1, 1.8, 1], vs=[0, 0, 1, 0.7])
p = plot(legend=:topright, framestyle=:box, grid=false, title="LinearLayeredModel example", xlabel="Radius (km)", ylabel="Velocity (km/s)")
for (f,s,l) in zip((vp, vs), (:vp, :vs), ("Vp", "Vs"))
    plot!(p, r->f(m, r), 0, surface_radius(m), label=l, line=2)
    scatter!(p, m.r, getfield(m, s), label="Radial nodes ($l)")
end
savefig(p, "LinearLayeredModel_example.svg")
nothing
```
![Example of a LinearLayeredModel](LinearLayeredModel_example.svg)

Use [`discontinuities`](@ref) to find the radii and layer indices of
the discontinuities in a `LinearLayeredModel`.

```@docs
LinearLayeredModel
```

See the [constructor](@ref SeisModels.LinearLayeredModel()) for details on creating
`LinearLayeredModel`s.

### `PREMPolyModel`

`PREMPolyModel`s are parameterised (like PREM) by polynomials in each layer,
considering the normalised radius $x=r/a$, where $a$ is the surface radius.
Clearly, `SteppedLayeredModel`s and `LinearLayeredModel`s are specific
cases of a polynomial parameterisation.

For example, PREM’s outer core density in g/cm³ is given by
$$12.5815 - 1.2638x - 3.6426x^2 - 5.5281x^3\,.$$

One could imagine a single-layer planet whose velocities are described by
```math
V_\mathrm{P}(x) = 3 - 2x^2 - x^3\,,
```
where $V_\mathrm{S} = V_\mathrm{P}/\alpha$, and $\alpha = 1.7$.

```@repl example
radii = [500]    # Layer tops in km
vps =   reshape( # Ensure we have a matrix with 1 column and 4 rows
        [3       # x^0 coefficients in km/s
         0       # x^1
        -2       # x^2
        -1],     # x^3
         :, 1)
α = 1.7
vss = vps./α

m = PREMPolyModel(r=radii, vp=vps, vs=vss)
```

Note that the polynomial coefficients appear as the row of a matrix,
and each layer occupies a column.

```@repl example
using Plots: plot, plot!, vline!, savefig
p = plot(xlabel="Radius (km)", ylabel="Velocity (km/s) or Density (g/cm^3)", title="PREMPolyModel example", framestyle=:box, grid=false)
for (f, name) in zip((vp, vs), ("Vp", "Vs"))
    plot!(p, r->f(m, r), 0, surface_radius(m), label=name)
end
vline!(p, m.r, line=(0.5,:dash), label="Layer tops")
savefig("PREMPolyModel_example.svg")
```
![Example of a PREMPolyModel](PREMPolyModel_example.svg)

```@docs
PREMPolyModel
```

See the [constructor](@ref SeisModels.PREMPolyModel()) for details on creating
`PREMPolyModel`s.

#### Attenuation
Attenuation in `PREMPolyModels` is specified as in PREM (see equation
(3) in section 6 on page 309).  That means, the
shear and bulk quality factors ($Q_\mu$ and $Q_\kappa$ respectively)
control the effective velocities, which will vary by wave frequency.

By default, velocities for `PREMPolyModel`s are returned at the reference
frequency, which can be found with [`reffrequency`](@ref) and is given
in Hz.  To obtain velocities at other frequencies, just pass a target
frequency to the relevant function via the `freq` keyword argument.

In this example, we ask for the Voigt average isotropic shear-wave
velocity at radius 4500 km for a frequency of 0.01 Hz (equally,
100 s period):
```@repl example
vs(PREM, 4500, freq=0.01)
```


## Conversion
Conversion can be performed between any of the above model types.  Note
that in some cases (e.g., `SteppedLayeredModel` to `LinearLayeredModel`
or `LinearLayeredModel` to `PREMPolyModel`) no information is lost, whilst
in others (e.g., `PREMPolyModel` to `LinearLayeredModel` or to
`SteppedLayeredModel`) the converted model can only be an approximation
to the original.

To convert between model types, simply pass one model to the constructor
of another.  For example:

```@repl example
m = PREM;
m2 = LinearLayeredModel(PREM);
typeof(m2)
≈(mass(m), mass(m2))
```

See the [conversion section](@ref conversion) for details.
