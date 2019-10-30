"""
# `MOON_WEBER_2011`

1-D layered model of the Moon, containing density.

## Reference

Weber, R.C.,, Lin, P.-Y.,, Garnero, E.J., Williams, Q., Lognonné, P., 2011.
Seismic Detection of the Lunar Core. Science, 331, 309–312.
https://doi.org/10.1126/science.1199375
"""
const MOON_WEBER_2011 = SteppedLayeredModel(1737.1, 10,
    # radius
    1737.1 .- reverse([  0.0,   1.0,   15.0,   40.0,  238.0,
                       488.0, 738.0, 1257.1, 1407.1, 1497.1]),
    # Vp
    reverse([1.0, 3.2, 5.5, 7.7, 7.8, 7.6, 8.5, 7.5, 4.1, 4.3]),
    # Vs
    reverse([0.5, 1.8, 3.2, 4.4, 4.4, 4.4, 4.5, 3.2, 0.0, 2.3]),
    # ρ
    reverse([2.6, 2.7, 2.8, 3.3, 3.4, 3.4, 3.4, 3.4, 5.1, 8.0]),
    # Anisotropy
    false, Float64[], Float64[], Float64[], Float64[], Float64[],
    # Attenuation
    false, Float64[], Float64[])
