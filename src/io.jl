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
    voigt_kappa_mu(A, C, L, N, F) -> κ, μ

Return the Voigt average bulk modulus `κ` and shear modulus `μ` given the Love parameters
of a transversely isotropic medium.  If A, C, etc., contain density, then κ and μ will be
true moduli, whilst if the Love parameters are provided as velocities, then
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
    