# Input and output routines

"Maximum polynomial degree of Mineos polynomial model files"
MINEOS_MAX_DEGREE = 4

"""
    read_mineos(file) -> m::LinearLayeredModel

Read a 1D seismic model from a Mineos tabular format file.

### References:
- https://geodynamics.org/cig/software/mineos/mineos-manual.pdf
"""
function read_mineos(io::IO)
    name = readline(io)
    ifanis, tref, ifdeck = parse.((Int, Float64, Int), split(readline(io)))
    aniso = ifanis == 1
    ifdeck != 1 && error("File '$file' is not in tabular format (`ifdeck` is `$ifdeck`)")
    n, nic, noc = parse.(Int, split(readline(io)))
    d = readdlm(io)
    r = d[:,1]./1e3
    a = maximum(r)
    rho = d[:,2]./1e3
    if aniso
        vpv = d[:,3]./1e3
        vsv = d[:,4]./1e3
        vph = d[:,7]./1e3
        vsh = d[:,8]./1e3
        eta = d[:,9]
        voigt_vels = voigt_velocities.(vpv, vsv, vph, vsh, eta)
        vp = [v[1] for v in voigt_vels]
        vs = [v[2] for v in voigt_vels]
    else
        vp = d[:,3]./1e3
        vs = d[:,4]./1e3
        vpv = vph = vsv = vsh = eta = []
    end
    Qκ = d[:,5]
    Qμ = d[:,6]
    attenuation = !(all(isequal(0), Qκ) && all(isequal(0), Qμ))
    if !attenuation
        Qμ = []
        Qκ = []
    end
    LinearLayeredModel(a, n, r, vp, vs, rho, aniso, vph, vpv, vsh, vsv, eta,
                       attenuation, Qμ, Qκ)
end

read_mineos(file::AbstractString) = (open(file, "r") do io; read_mineos(io); end)

# TODO: Remove when bumping minor or major version
Base.@deprecate write_mineos(model::LinearLayeredModel, file::AbstractString, args...) write_mineos(file, model, args...)

"""
    write_mineos(file, m::LinearLayeredModel, freq=1.0, title="Model from SeisModels.jl")
    write_mineos(io::IO, m, args...)

Save a `SeisModel1D` as a Mineos 'tabular' format file to `file` on disk or
to an `io`.  Optionally supply the reference
frequency `freq` in Hz (which defaults to 1 Hz) and a `title`.

Note that Mineos has two types of model file; one parameterised by radial knots
with values interpolated linearly between them ('tabular' format),
and one with PREM-style polynomials ('polynomial' format).
Mineos uses the tabular format internally even if a polynomial file is read in,
and testing suggests that polynomial files are not used correctly by
the [CIG](https://geodynamics.org/cig/software/mineos/) version of Mineos,
so it is not recommended to use them at all.  Consequently, SeisModels does
not support writing polynomial format files.

## Reference
- Mineos manual, https://geodynamics.org/cig/software/mineos/mineos-manual.pdf
"""
function write_mineos(io::IO, m::LinearLayeredModel, freq=1.0, title="Model from SeisModels.jl")
    _check_mineos_title(title)
    ifanis = isanisotropic(m) ? 1 : 0
    tref = freq
    ifdeck = 1
    N = _check_mineos_num_layers(m)
    nic, noc = core_interface_layers(m)
    rho, vp, vs = m.density.*1e3, m.vp.*1e3, m.vs.*1e3
    if m.aniso
        vpv, vph, vsv, vsh, eta = m.vpv.*1e3, m.vph.*1e3, m.vsv.*1e3, m.vsh.*1e3, m.eta
    else
        vpv = vph = vp
        vsh = vsv = vs
        eta = ones(m.n)
    end
    if hasattenuation(m)
        Qμ, Qκ = m.Qμ, m.Qκ
    else
        Qμ = Qκ = zeros(m.n)
    end
    println(io, title)
    println(io, ifanis, " ", tref, " ", ifdeck)
    println(io, N, " ", nic, " ", noc)
    writedlm(io, [m.r*1e3 rho vpv vsv Qκ Qμ vph vsh eta])
    nothing
end

write_mineos(file::AbstractString, m::LinearLayeredModel, args...) =
    (open(file, "w") do io; write_mineos(io, m, args...); end)

"""
    _check_mineos_title(title)

Issue a warning if `title` is longer than Mineos's maximum title length
(80 characters).
"""
_check_mineos_title(title) = length(title) > 80 &&
        @warn("Mineos model files can have titles only 80 characters long " *
              "('$title' is $(length(title)) characters)")

"""
    _check_mineos_num_layers(m::SeisModel) -> N

Issue a warning if the model `m` contains more than Mineos's maximum
number of layers (350), and return the number of layers `N`.
"""
_check_mineos_num_layers(m::SeisModel) = (N = m.n;
    N > 350 && @warn("Mineos model files are limited to N ≤ 350 (have $N layers)"); N)

"""
    _check_mineos_max_poly_degree(m::PREMPolyModel)

Throw an error if the model `m` has any fields with polynomial degree greater
than Mineos's maximum degree (4).
"""
function _check_mineos_max_poly_degree(m::PREMPolyModel)
    # Mineos only deals with up to degree-four polynomials with 5 coefficients
    max_degree = 4
    for field in (:vp, :vs, :density, :vph, :vpv, :vsh, :vsv, :eta, :Qμ, :Qκ)
        v = getfield(m, field)
        if ndims(v) > 1 && size(v, 1) > MINEOS_MAX_DEGREE + 1
            throw(ArgumentError("maximum degree polynomial in Mineos polynomial files is 4"))
        end
    end
    nothing
end

"""
    read_tvel(file) -> model::LinearLayeredModel
    read_tvel(io::IO) -> model::LinearLayeredModel

Read a `model` from a file in the ttimes 'tvel' format.  This format
specifies a `LinearLayeredModel` with only isotropic Vp, Vs and
density.

The two header lines are ignored but must be present.
"""
function read_tvel(io::IO; file=nothing)
    lines = readlines(io)
    length(lines) > 2 || error("input .tvel string has too few lines")
    r, Vp, Vs, rho = Float64[], Float64[], Float64[], Float64[]
    for i in 3:length(lines)
        tokens = split(lines[i])
        length(tokens) >= 4 || error("too few columns at line $i of " *
            (file === nothing ? "input" : "file '$file'"))
        rr, vvpp, vvss, ρρ = parse.(Float64, tokens)
        push!(r, rr)
        push!(Vp, vvpp)
        push!(Vs, vvss)
        push!(rho, ρρ)
    end
    R = r[end]
    r .= R .- r
    reverse!(r)
    reverse!(Vp)
    reverse!(Vs)
    reverse!(rho)
    LinearLayeredModel(r=r, vp=Vp, vs=Vs, density=rho)
end

read_tvel(file::AbstractString) = (open(file, "r") do io
    read_tvel(io, file=file); end)

"""
    write_tvel(file, model::LinearLayeredModel; comment1="SeisModels.jl - P", comment2="SeisModels.jl - S")
    write_tvel(io::IO, model; kwargs...)

Write a `model` to `file` or an arbitrary `IO` object
in the `ttimes` 'tvel' format.

Optionally supply the two comment lines, `comment1` and `comment2`
which form the header of the file.
"""
function write_tvel(io::IO, model::LinearLayeredModel; file=nothing,
    comment1="SeisModels.jl - P", comment2="SeisModels.jl - S")
    println(io, comment1)
    println(io, comment2)
    dep = surface_radius(model) .- reverse(model.r)
    writedlm(io, [dep reverse(model.vp) reverse(model.vs) reverse(model.density)])
end

write_tvel(file::AbstractString, model::LinearLayeredModel; kwargs...) =
    (open(file, "w") do io; write_tvel(io, model; file=file, kwargs...); end)

"""
    core_interface_layers(m) -> i_icb, i_cmb

Return the indices of the inner-core boundary, `i_icb`, and the core-mantle boundary,
`i_cmb`, for the model `m` based on the following assumptions:

1. The outer core is a liquid and has Vs = 0.
2. There is only one outer core
3. The inner core is solid and has Vs > 0.
4. There is only one inner core
5. The mantle above the outer core is solid and has Vs > 0.

These numbers refer to the layer which is the **top** of the inner/outer core.
"""
function core_interface_layers(m)
    i = 1
    i_icb = i_cmb = 0
    got_icb = false
    while i <= m.n
        if !got_icb && isapprox(m.vs[i], 0, atol=eps(eltype(m.vs)))
            got_icb = true
            i > 1 || error("Unexpectedly found Vs = 0 at centre of model")
            i_icb = i - 1
        end
        if got_icb && m.vs[i] > eps(eltype(m.vs))
            i_cmb = i - 1
            return i_icb, i_cmb
        end
        i += 1
    end
end

"""
    core_interface_layers(m::PREMPolyModel) -> i_icb, i_cmb

Return the ICB and CMB indices for a `PREMPolyModel`, using the following
heuristic:

1. The ICB and CMB are located at boundaries between layers.
2. The velocity 1 m above the inner core is liquid and has Vs = 0.
3. There is a solid inner core at the centre of the model.
4. The mantle is above a liquid, above a solid, and is itself solid.
"""
function core_interface_layers(m::PREMPolyModel)
    i = 1
    i_icb = i_cmb = 0
    got_icb = false
    while i <= m.n
        # Radii one meter either side of a layer boundary
        r_plus_one = m.r[i] + 1e-3
        if !got_icb && isapprox(vs(m, r_plus_one), 0, atol=eps(eltype(m.vs)))
            got_icb = true
            i_icb = i
        end
        if got_icb && vs(m, r_plus_one) > eps(eltype(m.vs))
            i_cmb = i
            return i_icb, i_cmb
        end
        i += 1
    end
end
