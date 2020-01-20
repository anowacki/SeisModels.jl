---
title: "SeisModels.jl: A Julia package for models of the Earth's interior"
tags:
  - Julia
  - geophysics
  - seismology
  - Earth
authors:
  - name: Andy Nowacki
    orcid: 0000-0001-7669-7383
    affiliation: 1
affiliations:
  - name: School of Earth and Environment, University of Leeds
    index: 1

date: 13 January 2020
bibliography: paper.bib
---



# Summary

The interior structure of the Earth and other quasi-spherical planets
reflects their history, development, present-day state and dynamics.  Hence
inferring their structure is of great interest, and requires the ability
to rapidly compute their bulk properties such as total mass, moment of inertia,
gravitational acceleration, interior pressure, and so on, in order to compare
model predictions with data.  ``SeisModels.jl`` is a Julia package designed to
allow such computation.

``SeisModels.jl`` can be used to represent arbitrary models of quasi-spherical
bodies.  In its current release, radially-symmetric bodies can be
represented easily using spherical shells whose properties are linearly varying with radius, constant with radius, or parameterised
by a set of polynomial coefficients.  Density,
elastic isotropic velocity (P-wave and S-wave), elastic radially-anisotropic
velocity and attenuation are supported as basic model parameters,
permitting the computation at arbitrary radius of gravity, pressure,
mass, and so on.  Inbuilt models include PREM [@dziewonski:1981],
AK135 [@kennett:1995] and iasp91 [@kennett:1991] for the Earth,
and the Moon model of @weber:2011.

Due to Julia's multiple dispatch programming paradigm, extending
``SeisModels.jl`` to represent bodies of different types or ones parameterised
differently to the inbuilt types is simple and incurs no performance
penalty.  The software has been benchmarked against independent implementations
and computes properties quickly enough to test millions of models against
data, enabling robust parameter searches or Monte Carlo sampling.
For the example of working with PREM on a desktop computer with an Intel
E5-1650 CPU, it takes only 12 ms to calculate the pressure at the centre
of the Earth, whilst surface gravity is returned within 5 µs.

``SeisModels.jl`` has been used for research projects and teaching,
and provides the ability to read and write files in appropriate formats
to calculate normal mode eigenfrequencies and -functions using the
[Mineos](https://geodynamics.org/cig/software/mineos/) software, and
ray-theoretical travel times using the Java TauP [@crotwell:1999]
and Python Obspy [@obspy] packages. ``SeisModels.jl`` is used by the
[``Mineos.jl``](https://github.com/anowacki/Mineos.jl)
Julia module to directly compute normal mode properties for planetary
models.

Documentation is available at a
[dedicated website](https://anowacki.github.io/SeisModels.jl/stable).
The package comes with a test set which can be easily run by
users—after adding the package, simply run
`import Pkg; Pkg.test("SeisModels")` in Julia.


# Acknowledgements

The author was supported by a NERC standard grant (NE/R001154/1).
The Deep Earth Research Group at the School of Earth and Environment, University of Leeds is acknowledged for inspiring the creation of this
package, and Andrew Walker for the publication of it.

# References

