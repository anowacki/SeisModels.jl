# Model input/output

SeisModels contains the ability to read and write [Mineos](https://geodynamics.org/cig/software/mineos/) files in tabular format.
These by definition are [`SteppedLayeredModel`](@ref)s.

Use the commands [`read_mineos`](@ref) and [`write_mineos`](@ref)
to respectively load and save `SteppedLayeredModel`s in this format.
