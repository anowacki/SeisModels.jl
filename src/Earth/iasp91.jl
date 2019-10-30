"""
# `IASP91`

1-D model of velocity.

**N.B.** iasp91 does not contain density.

This version of the model uses the original coefficients and is a
`PREMPolyModel`  as published.  Other implementations available
may use a layered piecewise linear model and/or add density from
other sources.

## Reference

Kennett  B.L.N.  Engdahl  E.  1991. Traveltimes for global earthquake location
and phase identification. Geophys J Int 105  429–465.
https://doi.org/10.1111/j.1365-246X.1991.tb06724.x
"""
const IASP91 = PREMPolyModel(6371.0, 11,
	# r
	[  1217.1,     3482,     3631,     5611,     5711,     5961,     6161,     6251,      6336, 6351, 6371],
	# Vp
	[11.24094  10.03904  14.49470  25.14860  25.96984  29.38896  30.78765  25.41389   8.787541   6.5   5.8;
	  0.0       3.75665  -1.47089 -41.15380 -16.93412 -21.40656 -23.25415 -17.69722  -0.749530   0.0   0.0;
	 -4.09689 -13.67046   0.0      51.99320   0.0       0.0       0.0       0.0       0.0        0.0   0.0;
	  0.0       0.0       0.0     -26.60830   0.0       0.0       0.0       0.0       0.0        0.0   0.0],
	# Vs
	[ 3.56454   0.0       8.16616  12.93030  20.76890  17.70732  15.24213   5.75020   6.706231   3.75  3.36;
	  0.0       0.0      -1.58206 -21.25900 -16.53147 -13.50652 -11.08552  -1.27420  -2.248585   0.0   0.0;
	 -3.45241   0.0       0.0      27.89880   0.0       0.0       0.0       0.0       0.0        0.0   0.0;
	  0.0       0.0       0.0     -14.10800   0.0       0.0       0.0       0.0       0.0        0.0   0.0],
	# rho
	[],
	# aniso
	false, [], [], [], [], [],
	# attenuation
	false, [], []
	)
