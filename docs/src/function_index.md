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
discontinuities
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
read_tvel
write_tvel
```

## [Conversion](@id conversion)
```@docs
LinearLayeredModel(::PREMPolyModel)
LinearLayeredModel(::SteppedLayeredModel)
PREMPolyModel(::LinearLayeredModel)
PREMPolyModel(::SteppedLayeredModel)
SteppedLayeredModel(::Union{LinearLayeredModel, PREMPolyModel})
```

## Construction
```@docs
LinearLayeredModel()
PREMPolyModel()
SteppedLayeredModel()
```