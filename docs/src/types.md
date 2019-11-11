# Model types

## Type hierarchy
The module defines the SeisModel type and subtypes of this specify the kind of
model (i.e., symmetry, nature of basis function, etc.).

The current type hierarchy is as follows:

- `SeisModel`: Abstract type of all seismic models
  - `SeisModel1D`: Abstract type of all one-dimensonal (radially-symmetric)
    seismic models
    - `SteppedLayeredModel`: Type defined by layers with contstat properties.
    - `LinearLayeredModel`: Type defined by layers with constant gradient.
    - `PREMPolyModel`: Type defined by polynomials.

As you can see, there are currently three types of models implemented, all
1D models, with polynomial, linear or constant basis within each layer.

### `SteppedLayeredModel`

A `SteppedLayeredModel` is composed of a series of layers with constant
properties in each layer.  The radii define the top of each layer.

An example of this is Weber et al.’s (2011) moon model:
```@eval
using SeisModels, Plots
p = plot(legend=:bottom, framestyle=:box, grid=false, title="SteppedLayeredModel example: the moon", xlabel="Radius (km)", ylabel="Velocity (km/s) or density (g/cm^3)")
for (f,l) in zip((vp, vs, SeisModels.density), ("Vp", "Vs", "Density"))
    plot!(p, r->f(MOON_WEBER_2011, r), 0, surface_radius(MOON_WEBER_2011), label=l, line=2)
end
vline!(p, MOON_WEBER_2011.r, label="Layer tops", l=(0.5, :dash))
savefig(p, "SteppedLayeredModel_example.svg")
```
![Example of a SteppedLayeredModel](SteppedLayeredModel_example.svg)

```@docs
SteppedLayeredModel
```

### `LinearLayeredModel`

A `LinearLayeredModel` is described by radial nodes, between which
the properties vary linearly.  Discontinutities are added by setting
two subsequent nodes to have the same radii but different properties.

For example, a two-layer model for a planetesimal with radius 1 km
and a liquid core at 0.3 km, with linearly varying properties, could
be given by:

```@repl example
using SeisModels
m = LinearLayeredModel(1.0, 4, [0, 0.3, 0.3, 1], [1, 1, 1.8, 1], [0, 0, 1, 0.7], [], false, [], [], [], [], [], false, [], [])
```

You can see that the nodes are used to define a constant lower layer
and linearly-varying upper layer:
```@eval
using Plots, SeisModels
m = LinearLayeredModel(1.0, 4, [0, 0.3, 0.3, 1], [1, 1, 1.8, 1], [0, 0, 1, 0.7], [], false, [], [], [], [], [], false, [], [])
p = plot(legend=:topright, framestyle=:box, grid=false, title="LinearLayeredModel example", xlabel="Radius (km)", ylabel="Velocity (km/s)")
for (f,s,l) in zip((vp, vs), (:vp, :vs), ("Vp", "Vs"))
    plot!(p, r->f(m, r), 0, surface_radius(m), label=l, line=2)
    scatter!(p, m.r, getfield(m, s), label="Radial nodes ($l)")
end
savefig(p, "LinearLayeredModel_example.svg")
nothing
```
![Example of a LinearLayeredModel](LinearLayeredModel_example.svg)

```docs
LinearLayeredModel
```

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
isaniso = hasQ = false
rhos = vphs = vpvs = vshs = vsvs = etas = Qmu = Qkappa = []

m = PREMPolyModel(500, 1, radii, vps, vss, rhos, isaniso, vphs, vpvs, vshs, vsvs, etas, hasQ, Qmu, Qkappa)
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
