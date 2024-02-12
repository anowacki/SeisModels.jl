# Generate the SC-REM models

let data_path = joinpath(@__DIR__, "..", "..", "data", "Earth", "SC_REM")
    for percent in (25, 50, 75), high_low in ("high", "low")
        model_file = joinpath(data_path, "screm-$(percent)p-$(high_low).dat")
        model_data = readdlm(model_file, ',', Float64; skipstart=1)

        r = model_data[:,1]
        ρ = model_data[:,3]
        vp = model_data[:,4]
        vs = model_data[:,5]
        Qκ = model_data[:,6]
        Qμ = model_data[:,7]

        lower_upper = high_low == "high" ? "Upper" : "Lower"

        model = Symbol(:SC_REM_, percent, :_, uppercase(high_low))
        model_name = String(model)

        @eval begin
            """
            # `$($model_name)`


            !!! warning "Retracted model"
                The $($model_name) model has been retracted by the authors and
                is only retained in SeisModels to retain compatibility.
                Future major versions of SeisModels may remove this model.

            $($lower_upper) bound of the $($percent)% credible interval of the SC-REM model of
            Kemper et al. (2023).

            See also:
            [`SC_REM_25_LOW`](@ref),
            [`SC_REM_25_HIGH`](@ref),
            [`SC_REM_50_LOW`](@ref),
            [`SC_REM_50_HIGH`](@ref),
            [`SC_REM_75_LOW`](@ref),
            [`SC_REM_75_HIGH`](@ref).

            ## References

            - Retraction of: Self-consistent models of Earth's mantle and core from long-period seismic and tidal constraints,
              Geophysical Journal International, Volume 236, Issue 3, March 2024, Page 1439,
              https://doi.org/10.1093/gji/ggae004
            - Kemper, J., Khan, A., Helffrich, G., van Driel, M., Giardini, D., 2023.
              Self-consistent models of Earth’s mantle and core from long-period seismic and tidal constraints.
              Geophys J Int 235, 690–717.
              https://doi.org/10.1093/gji/ggad254
            - Kemper, J., Khan, A., Helffrich, G., Giardini, D., 2023.
              Self-consistent models of Earth's mantle and core from long-period seismic and tidal constraints (0.2) [Data set].
              Zenodo.
              https://doi.org/10.5281/zenodo.10037465
            """
            const $model = LinearLayeredModel(
                r = $r,
                vp = $vp,
                vs = $vs,
                density = $ρ,
                Qκ = $Qκ,
                Qμ = $Qμ
            )
        end
    end
end
