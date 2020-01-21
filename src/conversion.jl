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
    Vp, Vs = vp.(m, r_eval), vs.(m, r_eval)
    Rho = hasdensity(m) ? density.(m, r_eval) : []
    Vpv, Vph, Vsv, Vsh, Eta = isanisotropic(m) ?
        (vpv.(m, r_eval), vph.(m, r_eval), vsv.(m, r_eval), vsh.(m, r_eval), eta.(m, r_eval)) :
        ([], [], [], [], [])
    qκ, qμ = hasattenuation(m) ? (Qκ.(m, r_eval), Qμ.(m, r_eval)) : ([], [])
    LinearLayeredModel(R, length(radii), radii, Vp, Vs, Rho,
        isanisotropic(m), Vph, Vpv, Vsh, Vsv, Eta,
        hasattenuation(m), qμ, qκ)
end

"Set of fields we copy to the LinearLayeredModel struct"
const LinearLayeredModel_CONVERSION_FIELDS =
    (x for x in fieldnames(LinearLayeredModel)
     if x ∉ (:a, :n, :r, :aniso, :attenuation))

"""
    LinearLayeredModel(m::SteppedLayeredModel)

Convert the model `m` into a `LinearLayeredModel`, with radial
knots spaced at the bottom and top of each layer in `m`.
"""
function LinearLayeredModel(m::SteppedLayeredModel)
    R = surface_radius(m)
    d = Dict{Symbol,Vector{Float64}}()
    radii = Iterators.flatten((i == 1 ? 0.0 : m.r[i-1], m.r[i]) for i in eachindex(m.r))
    for field in LinearLayeredModel_CONVERSION_FIELDS
        val = getfield(m, field)
        isempty(val) && continue
        d[field] = repeat(val, inner=2)
    end
    d[:r] = collect(radii)
    LinearLayeredModel(; d...)
end

"Set of fields we copy to the PREMPolyModel struct"
const PREMPolyModel_CONVERSION_FIELDS = (x for x in fieldnames(PREMPolyModel)
                                         if x ∉ (:r, :a, :n, :fref, :aniso, :attenuation))

"""
    PREMPolyModel(m::SteppedLayeredModel, order=0; fref=1.0)

Convert the model `m` into a `PREMPolyModel`, with maximum polynomial
`order`.  The reference frequency `fref` in Hz can be given.
"""
function PREMPolyModel(m::SteppedLayeredModel, order=0; fref=1.0)
    order < 0 &&
        throw(ArgumentError("cannot construct polynomials with negative order"))
    nlayers = length(m.r)
    d = Dict{Symbol,Any}()
    for field in PREMPolyModel_CONVERSION_FIELDS
        val = getfield(m, field)
        isempty(val) && continue
        if val isa AbstractArray && order > 0
            arr = fill(0.0, order + 1, nlayers)
            arr[1,:] .= val
            d[field] = arr
        else
            d[field] = val
        end
    end
    d[:fref] = fref
    d[:r] = m.r
    PREMPolyModel(; d...)
end

"""
    PREMPolyModel(m::LinearLayeredModel, order=1; fref=1.0)

Convert the model `m` into a `PREMPolyModel`, with maximum polynomial
`order`.  The reference frequency `fref` in Hz can be given.
"""
function PREMPolyModel(m::LinearLayeredModel, order=1; fref=1.0)
    order < 1 &&
        throw(ArgumentError("minimum polynomial order is 1 when " *
                            "converting `LinearLayeredModel`s"))
    d = Dict{Symbol,Any}()
    R = surface_radius(m)
    # Indices of nodes on top of finite layers (i.e., not discontinuity layers)
    layer_inds = [i for i in 2:length(m.r) if m.r[i] > m.r[i-1]]
    nlayers = length(layer_inds)
    for field in PREMPolyModel_CONVERSION_FIELDS
        val = getfield(m, field)
        isempty(val) && continue
        if val isa AbstractArray
            d[field] = fill(0.0, order + 1, nlayers)
            for (i, ind) in enumerate(layer_inds)
                d[field][1:2,i] .= _fit_prem_polynomial(R,
                    m.r[ind-1], m.r[ind], val[ind-1], val[ind])
            end
        else
            d[field] = val
        end
    end
    d[:fref] = fref
    d[:r] = m.r[layer_inds]
    PREMPolyModel(; d...)
end

"""
    _fit_prem_polynomial(a, r1, r2, v1, v2) -> (a0, a1)

Return PREM polynomial coefficients (where `a0` is the constant and `a1`
the linear coefficient) for the line connecting two values `v1` and
`v2` at radii `r1` and `r2`.  The radii are unnormalised compared to
the planetary radius `a`.
"""
function _fit_prem_polynomial(a, r1, r2, v1, v2)
    x1, x2 = r1/a, r2/a
    a1 = (v2 - v1)/(x2 - x1)
    a0 = v1 - a1*x1
    a0, a1
end

"Set of fields we copy into the SteppedLayeredModel struct"
const SteppedLayeredModel_CONVERSION_FIELDS =
    (x for x in fieldnames(SteppedLayeredModel)
     if x ∉ (:r, :a, :n, :fref, :aniso, :attenuation))

"""
    SteppedLayeredModel(m::Union{LinearLayeredModel, PREMPolyModel}, spacing=10)

Convert the model `m` into a `SteppedLayeredModel`, with layer tops
spaced a maximum of `spacing` km apart.  Where layer tops of `m` lie
between the regular set of `spacing` layer tops, extra layers are
inserted.  Thus layer tops in the new model always coincide with
discontinuities in `m`.

Layers take the value of the properties at the midpoint of the
`spacing` km interval in the original `PREMPolyModel`.
"""
function SteppedLayeredModel(m::Union{LinearLayeredModel, PREMPolyModel}, spacing=10)
    spacing <= 0 && throw(ArgumentError("layer spacing must be positive"))
    d = Dict{Symbol,Any}()
    R = surface_radius(m)
    # Regular set of radii spaced `spacing` km apart
    spaced_radii = reverse(range(R, stop=spacing, step=-spacing))
    # New model layer top radii with m's layers interspersed.
    radii = unique!(sort!(vcat(m.r, spaced_radii)))
    # Radii from which to take layer properties in the new model
    eval_radii = [((i == 1 ? 0.0 : radii[i-1]) + radii[i])/2 for i in eachindex(radii)]
    for field in SteppedLayeredModel_CONVERSION_FIELDS
        val = getfield(m, field)
        isempty(val) && continue
        d[field] = evaluate.(m, field, eval_radii)
    end
    d[:r] = radii
    SteppedLayeredModel(; d...)
end
