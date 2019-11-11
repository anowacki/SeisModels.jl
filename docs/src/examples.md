# Examples

## Preamble

These commands all assume that you have done
```@repl example
using SeisModels
```
first.

## Basic properties

To compute the P-wave velocity at the base of the Earth’s mantle
(2750 km depth) in the iasp91, ak135 and PREM models, and compare
them all to PREM, you could do
```@repl example
vels = vp.((IASP91, AK135, PREM), 2750, depth=true)
100 .* (vels .- vels[end]) ./ vels[end]
```

(So they all agree to less than a quarter of a percent!)

To retrieve a density profile across the inner core boundary in
PREM:
```@repl example
r = 1200:10:1250;
density.(PREM, r)
```

To get the planetary radius in km of the moon used in Weber et al.’s
(2011) model:
```@repl example
surface_radius(MOON_WEBER_2011)
```

The shear modulus in the transition zone:
```@repl example
shear_modulus(AK135, 600, depth=true)
```

## Derived properties

If a model includes density, then things like gravity, moment of inertia,
and so on, can be easily calculated.

PREM used as constraining data a couple of gross Earth values:
the mass should be $5.974\times10^{24}$kg, and the ratio
$I/MR^2$ should be 0.3308.  Let’s check how Dziewoński and Anderson
did:
```@repl example
≈(mass(PREM, 0, depth=true), 5.974e24, atol=0.001e24)
≈(moment_of_inertia(PREM)/(mass(PREM, 0, depth=true)*(1e3*surface_radius(PREM))^2), 0.3308, atol=0.0001)
```

To calculate a lookup table of pressure in ak135 for use with another
program, you could do:
```@repl example
using DelimitedFiles: writedlm
r = 0:10:6371
p = pressure.(AK135, r)
writedlm("lookup_table.txt", [r p])
```
