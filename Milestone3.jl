module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM()
xrm = equality_rm(co)
xmΛ = equality_mΛ(co)
xrec = time_recombination(co)
#co = ΛCDM(h=0.7, Neff=0, Ωb0=0.05, Ωc0=0.45, Yp=0, z_reion_H=NaN) # to compare with Hans
ks = [1e-0, 1e-1, 1e-2, 1e-3] / Mpc
#println("Tight coupling end (k = $(k*Mpc)/Mpc): ", format_time_variations(co, xtce))

# use spline points for plotting,
# but add nextra points between each of them
# TODO: make accessible to plotter?
# take an array of x values (e.g. spline points),
# then add nextra points between each of them
function extendx(x::Vector{Float64}, nextra::Integer)
    dx = diff(x)
    return sort(vcat(x, (x[1:end-1] .+ i/(nextra+1)*dx for i in 0:nextra)...))
end

if false
    series = [
        ("plots/Thetal0.pdf", Dict(:ylabel => L"\Theta_0"), [((x,y) -> y[Cosmology.i_Θl(0)](x), Dict(:linestyle => :solid, :label => nothing))]),
        ("plots/Thetal1.pdf", Dict(:ylabel => L"\Theta_1"), [((x,y) -> y[Cosmology.i_Θl(1)](x), Dict(:linestyle => :solid, :label => nothing))]),
        ("plots/Thetal2.pdf", Dict(:ylabel => L"\Theta_2"), [((x,y) -> y[Cosmology.i_Θl(2)](x), Dict(:linestyle => :solid, :label => nothing))]),

        ("plots/potentials.pdf", Dict(:ylabel => L"\{\Phi,\Psi\}"), [((x,y) -> y[Cosmology.i_Φ](x), Dict(:linestyle => :dash, :label => L"\Phi")), ((x,y) -> y[Cosmology.i_Ψ](x), Dict(:linestyle => :dot, :label => L"\Psi")), ((x,y) -> y[Cosmology.i_Φ](x) + y[Cosmology.i_Ψ](x), Dict(:linestyle => :solid, :label => L"\Phi+\Psi"))]),

        ("plots/velocity1.pdf",    Dict(:ylabel => L"\log_{10}|v|"),      [((x,y) -> log10(abs(   y[Cosmology.i_vc   ](x))), Dict(:linestyle => :solid, :label => L"v=v_c")),                          ((x,y) -> log10(abs(   y[Cosmology.i_vb   ](x))), Dict(:linestyle => :dash, :label => L"v=v_b"                           ))]),
        ("plots/velocity2.pdf",    Dict(:ylabel => L"\log_{10}|v|"),      [((x,y) -> log10(abs(-3*y[Cosmology.i_Θl(1)](x))), Dict(:linestyle => :solid, :label => L"v=v_\gamma=-3\Theta_1")),          ((x,y) -> log10(abs(-3*y[Cosmology.i_Nl(1)](x))), Dict(:linestyle => :dash, :label => L"v=v_\nu=-3\mathcal{N}_1"         ))]),
        ("plots/overdensity1.pdf", Dict(:ylabel => L"\log_{10}|\delta|"), [((x,y) -> log10(abs(   y[Cosmology.i_δc   ](x))), Dict(:linestyle => :solid, :label => L"\delta=\delta_c")),                ((x,y) -> log10(abs(   y[Cosmology.i_δb   ](x))), Dict(:linestyle => :dash, :label => L"\delta=\delta_b"                 ))]),
        ("plots/overdensity2.pdf", Dict(:ylabel => L"\log_{10}|\delta|"), [((x,y) -> log10(abs(+4*y[Cosmology.i_Θl(0)](x))), Dict(:linestyle => :solid, :label => L"\delta=\delta_\gamma=4\Theta_0")), ((x,y) -> log10(abs(+4*y[Cosmology.i_Nl(0)](x))), Dict(:linestyle => :dash, :label => L"\delta=\delta_\nu=4\mathcal{N}_0"))]),

        ("plots/Nl0.pdf", Dict(:ylabel => L"\mathcal{N}_0"), [((x,y) -> y[Cosmology.i_Nl(0)](x), Dict(:linestyle => :solid, :label => nothing))]),
        ("plots/Nl1.pdf", Dict(:ylabel => L"\mathcal{N}_1"), [((x,y) -> y[Cosmology.i_Nl(1)](x), Dict(:linestyle => :solid, :label => nothing))]),
        ("plots/Nl2.pdf", Dict(:ylabel => L"\mathcal{N}_2"), [((x,y) -> y[Cosmology.i_Nl(2)](x), Dict(:linestyle => :solid, :label => nothing))]),

        ("plots/ThetaPl0.pdf", Dict(:ylabel => L"\Theta^P_0"), [((x,y) -> y[Cosmology.i_ΘPl(0)](x), Dict(:linestyle => :solid, :label => nothing))]),
        ("plots/ThetaPl1.pdf", Dict(:ylabel => L"\Theta^P_1"), [((x,y) -> y[Cosmology.i_ΘPl(1)](x), Dict(:linestyle => :solid, :label => nothing))]),
        ("plots/ThetaPl2.pdf", Dict(:ylabel => L"\Theta^P_2"), [((x,y) -> y[Cosmology.i_ΘPl(2)](x), Dict(:linestyle => :solid, :label => nothing))]),
    ]

    # pre-compute callable splines once and for all (index and call as y1s[i_k][i_qty](x))
    y1s = [[x -> spline(x, k) for spline in Cosmology.perturbations_splines(co; tight=false)] for k in ks] # 2D (x, k) splines for full system
    y2s = [Cosmology.perturbations_mode(co, k; tight=false) for k in ks] # 1D (x) splines for each "exact" k for       full system

    for (filename, plotsettings, func_linesettings) in series
        println("Plotting $filename")
        plot(xlabel=L"x = \log a", xlims=(-15, 0), xticks=-25:1:5, legend_position=:topleft; plotsettings...)
        for (i, k) in enumerate(ks)
            xhor = time_horizon_entry(co, k)
            x = extendx(Cosmology.splinex(y2s[i][1]), 7) # plot with extra points between every spline point for more smoothness
            for (func, linesettings) in func_linesettings
                plot!(x, func.(x, Ref(y1s[i])), alpha=0.5, linewidth=1.0, color=i; linesettings..., label=nothing)
                plot!(x, func.(x, Ref(y2s[i])), alpha=1.0, linewidth=0.5, color=i; linesettings..., label=nothing)
                plot!([xhor], [func(xhor, y2s[i])], color=i, markerstrokecolor=i, markersize=1.0, markershape=:circle, label=nothing)
            end
            vline!([-21], color=i, label=L"k=10^{%$(Int(round(log10(k*Mpc))))}/\textrm{Mpc}") # label each k-value once
            #vline!([time_tight_coupling(co, k)], color=:gray, linestyle=:dash; linewidth=0.5, label=nothing)
            vline!([xrm, xmΛ], color=:gray, linestyle=:dash; linewidth=0.5, label=nothing)
            vline!([xrec], color=:gray, linestyle=:dot; linewidth=0.5, label=nothing)
        end

        # add quantity (if more than one so ambiguous) to legend with a black dummy plot
        for (func, linesettings) in func_linesettings
            vline!([-21], color=:black; linesettings...) # dummy outside plot area
        end
        savefig(filename)
    end
end

if false
    # plot Θ0 zoomed in (to check whether we capture the most rapid oscillations)
    filename = "plots/Thetal0_zoom.pdf"
    println("Plotting $filename")
    k = ks[2]
    y = Cosmology.perturbations_mode(co, k; tight=false)
    x = Cosmology.splinex(y[1])
    x = x[x .> -1.1]
    x = extendx(x, 20)
    func(x, y) = y[Cosmology.i_Θl(0)](x)
    plot(xlabel=L"x = \log a", ylabel=L"\Theta_0", xlims=(-1.0, 0), xticks=-5.0:1:0, legend_position=:topright)
    y1 = Cosmology.perturbations_mode(co, k; tight=true ) # 1D (x) splines for each "exact" k for tight+full system
    y2 = Cosmology.perturbations_mode(co, k; tight=false) # 1D (x) splines for each "exact" k for       full system
    plot!(x, func.(x, Ref(y1)), label=L"\texttt{ Tsit5}    \textrm{ method with tight coupling}")
    plot!(x, func.(x, Ref(y2)), label=L"\texttt{ KenCarp4} \textrm{ method without tight coupling}")
    savefig(filename)
end

# TODO: compare RK methods. plot max(yM1 - yM2) for each x and some k
if true
    ks = [1e-1, 1e-2, 1e-3] / Mpc
    filename = "plots/perturbation_methods.pdf"
    println("Plotting $filename")
    plot(xlabel=L"x = \log a", ylabel=L"\log_{10} \Big( \max_i{|\hat{y}^{\tiny{\textrm{KC4}}}_i-\hat{y}^{\tiny{\textrm{T5}}}_i|} \Big)", xlims=(-15, 0), xticks=-25:5:5, legend_position=:topleft)
    for (i_k, k) in enumerate(ks)
        y1s  = Cosmology.perturbations_mode(co, k; tight=true )
        y2s  = Cosmology.perturbations_mode(co, k; tight=false)
        x = extendx(Cosmology.splinex(y2s[1]), 20) # plot with extra points between every spline point for more smoothness
        y1s = [y1.(x) / maximum(abs.(y1.(x))) for y1 in y1s] # [i_q, i_x]
        y2s = [y2.(x) / maximum(abs.(y2.(x))) for y2 in y2s] # [i_q, i_x]
        dy = [abs.(y1 .- y2) for (y1, y2) in zip(y1s, y2s)]
        dy = [maximum([dy[i][j] for i in 1:length(dy)]) for j in 1:length(x)]
        plot!(x, log10.(dy); color=1+i_k, label=L"k=10^{%$(Int(round(log10(k*Mpc))))}/\textrm{Mpc}")
    end
    savefig(filename)
end

end
