# Model input/output

## Mineos files
SeisModels contains the ability to read and write [Mineos](https://geodynamics.org/cig/software/mineos/) files in tabular format.
These by definition are [`LinearLayeredModel`](@ref)s.

Use the commands [`read_mineos`](@ref) and [`write_mineos`](@ref)
to respectively load and save `LinearLayeredModel`s in this format.

[`PREMPolyModel`](@ref)s can be converted easily to `LinearLayeredModel`s:
see the section on [conversion](@ref PREMPolyModel) for details.

## `ttimes`-compatible `.tvel` files
You can also read and write '.tvel'-format files, as used by the
IASPEI [`ttimes`](http://rses.anu.edu.au/seismology/soft/ttsoft.html)
program, using [`read_tvel`](@ref) and [`write_tvel`](@ref).

'.tvel' files are also `LinearLayeredModel`s by definition—again,
other model types must be converted to this form before writing.
