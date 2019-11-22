# Conversion between model types

"""
    LinearLayeredModel(m::PREMPolyModel, spacing=20)

Convert the model `m` into a `LinearLayeredModel`, with radial
knots spaced a maximum of `spacing` km apart.
"""
function LinearLayeredModel(m::PREMPolyModel, spacing=20)
    R = surface_radius(m)
    # Radii for evaluation of model points at knots
    r_eval = Float64[]
    # Radii of knots, actually saved in new model
    radii = Float64[]
    for (ilayer, rtop) in enumerate(m.r)
        rtop_minus = prevfloat(rtop)
        # Bottom of layer
        rbot = ilayer == 1 ? 0.0 : m.r[ilayer - 1]
        rbot_plus = ilayer == 1 ? rbot : nextfloat(rbot)
        push!(radii, rbot)
        push!(r_eval, rbot_plus)
        r = rbot + spacing
        # Add knots
        while r < rtop
            push!(r_eval, r)
            push!(radii, r)
            r += spacing
        end
        # Top of layer
        push!(radii, rtop)
        push!(r_eval, rtop == R ? rtop : rtop_minus)
    end
    Vp, Vs, Rho = vp.(m, r_eval), vs.(m, r_eval), density.(m, r_eval)
    Vpv, Vph, Vsv, Vsh, Eta = isanisotropic(m) ?
        (vpv.(m, r_eval), vph.(m, r_eval), vsv.(m, r_eval), vsh.(m, r_eval), eta.(m, r_eval)) :
        ([], [], [], [], [])
    qκ, qμ = hasattenuation(m) ? (Qκ.(m, r_eval), Qμ.(m, r_eval)) : ([], [])
    LinearLayeredModel(R, length(radii), radii, Vp, Vs, Rho,
        isanisotropic(m), Vph, Vpv, Vsh, Vsv, Eta,
        hasattenuation(m), qμ, qκ)
end
