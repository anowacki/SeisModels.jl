# Index of functions

## Basic physical properties
```@docs
evaluate
vp
vs
density
vph
vpv
vsh
vsv
eta
Qμ
Qmu
Qκ
Qkappa
```

## Model enquiry
```@docs
depth
hasattenuation
hasdensity
hasreffrequency
isanisotropic
radius
reffrequency
surface_radius
```

## Derived model properties
```@docs
bulk_modulus
gravity
mass
moment_of_inertia
poissons_ratio
pressure
shear_modulus
surface_mass
youngs_modulus
```

## Input/output
```@docs
read_mineos
write_mineos
```

## Conversion
```@docs
LinearLayeredModel(::PREMPolyModel)
```

## Construction
```@docs
LinearLayeredModel()
PREMPolyModel()
SteppedLayeredModel()
```