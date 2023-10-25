# Test the implementation of SC_REM_* models by repdroducing Figure 4
# of Kemper et al. (2023)

using SeisModels
import Plots

function ribbon(xs, y1s, y2s)
    [xs; reverse(xs); xs[begin]], [y1s; reverse(y2s); y1s[begin]]
end

function ribbon_plot!(p, xs, y1s, y2s; kwargs...)
    x, y = ribbon(xs, y1s, y2s)
    Plots.plot!(p, Plots.Shape(y, x); kwargs...)
end

function ribbon_plot!(p, f, m1::SeisModels.SeisModel, m2::SeisModels.SeisModel, rs; kwargs...)
    vals1 = f.(m1, rs)
    vals2 = f.(m2, rs)
    ribbon_plot!(p, r, vals1, vals2; kwargs...)
end

ribbon_plot(args...; kwargs...) = ribbon_plot!(Plots.plot(), args...; kwargs...)

rs = 0:6371
prem_colour = :red

ps = Plots.Plot[]
for f in (density, vp, vs, Qμ)
    p = f == density ? Plots.plot() :
        f == Qμ ? Plots.plot(xaxis=:log10, xlim=(80, 1e5), legend=false) :
        Plots.plot(yticks=nothing, legend=false)
    push!(ps, p)
    for (low, high, color) in ((SC_REM_75_LOW, SC_REM_75_HIGH, :lightblue),
                               (SC_REM_50_LOW, SC_REM_50_HIGH, :blue),
                               (SC_REM_25_LOW, SC_REM_25_HIGH, :darkblue))
        ribbon_plot!(ps[end], f, low, high, rs;
            label="", fillcolor=color, title=f)
    end
    # Prem
    Plots.plot!(ps[end], f.(PREM, rs), rs, lc=prem_colour, ls=:dash, label="PREM")
end
Plots.plot(ps...; layout=(1, length(ps)), size=(800,400)) |> display
