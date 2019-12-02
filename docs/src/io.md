# Model input/output

SeisModels contains the ability to read and write [Mineos](https://geodynamics.org/cig/software/mineos/) files in tabular format.
These by definition are [`LinearLayeredModel`](@ref)s.

Use the commands [`read_mineos`](@ref) and [`write_mineos`](@ref)
to respectively load and save `LinearLayeredModel`s in this format.

[`PREMPolyModel`](@ref)s can be converted easily to `LinearLayeredModel`s:
see the section on [conversion](@ref PREMPolyModel) for details.
