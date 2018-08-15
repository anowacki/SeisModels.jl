# Input and output routines

"""
    read_mineos(file) -> ::SteppedLayeredModel

Read a `SteppedLayeredModel` in Mineos tabular format, which contains:

Line 1: Title
Line 2: anisotropy flag (0=false, 1=true), referencefrequency (Hz), tabular format flag (1=true)
Line 3: number of layers, inner core boundary index, core mantle boundary index
Lines 4+: radius (m), density (kg/m³), Vpv (m/s), Vsv (m/s), Qκ, Qμ, Vph (m/s),
          Vsh (m/s), η

### References:

1. https://geodynamics.org/cig/software/mineos/mineos-manual.pdf
"""
function read_mineos(file)
    open(file, "r") do f
        name = readline(f)
        ifanis, tref, ifdeck = parse.((Int, Float64, Int), split(readline(f)))
        aniso = ifanis == 1
        ifdeck != 1 && error("File '$file' is not in tabular format (`ifdeck` is `$ifdeck`)")
        n, nic, noc = parse.(Int, split(readline(f)))
        d = readdlm(f)
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
        SteppedLayeredModel(a, n, r, vp, vs, rho, aniso, vph, vpv, vsh, vsv, eta,
                            attenuation, Qμ, Qκ)
    end
end

"""
    save_mineos(m::SteppedLayeredModel, file, freq=1.0, title="Model from EarthModels.jl")

Save an `EarthModel1D` as a tabular Mineos-format file.  Supply the reference
frequency `freq` in Hz and a `title`.
"""
function write_mineos(m::SteppedLayeredModel, file, freq=1.0, title="Model from EarthModels.jl")
    length(title) > 80 &&
        warn("Mineos model files can have titles only 80 characters long " *
             "('$title' is $(length(title)) characters)")
    ifanis = m.aniso ? 1 : 0
    tref = freq
    ifdeck = 1
    N = m.n
    N > 350 && warn("Mineos model files are limited to N ≤ 350 (have $N layers)")
    nic, noc = core_interface_layers(m)
    rho, vp, vs = m.rho.*1e3, m.vp.*1e3, m.vs.*1e3
    if m.aniso
        vpv, vph, vsv, vsh, eta = m.vpv.*1e3, m.vph.*1e3, m.vsv.*1e3, m.vsh.*1e3, m.eta
    else
        vpv = vph = vp
        vsh = vsv = vs
        eta = ones(m.n)
    end
    if m.attenuation
        Qμ, Qκ = m.Qμ, m.Qκ
    else
        Qμ = Qκ = zeros(m.n)
    end
    open(file, "w") do f
        println(f, title)
        println(f, ifanis, " ", tref, " ", ifdeck)
        println(f, N, " ", nic, " ", noc)
        writedlm(f, [m.r*1e3 rho vpv vsv Qκ Qμ vph vsh eta])
    end
    nothing
end

"""
    core_interface_layers(m::SteppedLayeredModel) -> i_icb, i_cmb

Return the indices of the inner-core boundary, `i_icb`, and the core-mantle boundary,
`i_cmb`, for the model `m` based on the following assumptions:

1. The outer core is a liquid and has Vs = 0.
2. There is only one outer core
3. The inner core is solid and has Vs > 0.
4. There is only one inner core
5. The mantle above the outer core is solid and has Vs > 0.

These numbers refer to the layer which is the **top** of the inner/outer core.
"""
function core_interface_layers(m::SteppedLayeredModel)
    i = 1
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
    voigt_kappa_mu(A, C, L, N, F) -> κ, μ

Return the Voigt average bulk modulus `κ` and shear modulus `μ` given the Love parameters
of a transversely isotropic medium.  If A, C, etc., contain density, then κ and μ will be
true moduli, whilst if the Love parameters are provided as velocities, then these
are density-normalised 'moduli' with dimensions of L²⋅T².
"""
voigt_kappa_mu(A, C, L, N, F) = (4A + C + 4F - 4N)/9, (A + C - 2F + 5N + 6L)/15

"""
    voigt_velocities(vpv, vsv, vph, vsh, η) -> vp, vs

Return the Voigt average `vp` and `vs`, in the same units as the input, given four
velocities and the radial anisotropy parameter `η`.
"""
function voigt_velocities(vpv, vsv, vph, vsh, η)
    A = vph^2
    C = vpv^2
    L = vsv^2
    N = vsh^2
    F = η*(A - 2L)
    κ, μ = voigt_kappa_mu(A, C, L, N, F) # In km²/s²
    vp = sqrt(κ + 4/3*μ)
    vs = sqrt(μ)
    vp, vs
end
    