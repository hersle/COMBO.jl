module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM()
#co = ΛCDM(h=0.7, Neff=0, Ωb0=0.05, Ωc0=0.45, Yp=0, z_reion_H=NaN) # to compare with Hans
ks = [1e-3, 1e-2, 1e-1] / Mpc
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

if true
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
    y2s = [Cosmology.perturbations_mode(co, k; tight=true ) for k in ks] # 1D (x) splines for each "exact" k for tight+full system
    y3s = [Cosmology.perturbations_mode(co, k; tight=false) for k in ks] # 1D (x) splines for each "exact" k for       full system

    # plot Θ0 zoomed in (to check whether we capture the most rapid oscillations)
    filename = "plots/Thetal0_zoom.pdf"
    println("Plotting $filename")
    i = 3
    k = ks[i]
    y = Cosmology.perturbations_mode(co, k; tight=false)
    x = Cosmology.splinex(y[1])
    x = x[x .> -2.0]
    x = extendx(x, 20)
    func(x, y) = y[Cosmology.i_Θl(0)](x)
    plot(xlabel=L"x = \log a", ylabel=L"\Theta_0", xlims=(-1.0, 0), xticks=-2.0:0.1:0, legend_position=:topright)
    plot!(x, func.(x, Ref(y1s[i])), alpha=0.2, linewidth=1.5, color=i, label=nothing)
    plot!(x, func.(x, Ref(y2s[i])), alpha=0.5, linewidth=1.0, color=i, label=nothing)
    plot!(x, func.(x, Ref(y3s[i])), alpha=1.0, linewidth=0.5, color=i, label=L"k=%$(k*Mpc)/\textrm{Mpc}")
    savefig(filename)

    return

    for (filename, plotsettings, func_linesettings) in series
        println("Plotting $filename")
        plot(xlabel=L"x = \log a", xlims=(-20, 0), xticks=-25:5:5, legend_position=:topleft; plotsettings...)
        for (i, k) in enumerate(ks)
            x = extendx(Cosmology.splinex(y3s[i][1]), 7) # plot with extra points between every spline point for more smoothness
            for (func, linesettings) in func_linesettings
                plot!(x, func.(x, Ref(y1s[i])), alpha=0.2, linewidth=1.5, color=i; linesettings..., label=nothing)
                plot!(x, func.(x, Ref(y2s[i])), alpha=0.5, linewidth=1.0, color=i; linesettings..., label=nothing)
                plot!(x, func.(x, Ref(y3s[i])), alpha=1.0, linewidth=0.5, color=i; linesettings..., label=nothing)
            end
            vline!([-21], color=i, label=L"k = %$(k*Mpc) / \textrm{Mpc}") # label each k-value once
            vline!([time_tight_coupling(co, k)], color=:gray, linestyle=:dash; linewidth=0.5, label=nothing)
        end

        # add quantity (if more than one so ambiguous) to legend with a black dummy plot
        for (func, linesettings) in func_linesettings
            vline!([-21], color=:black; linesettings...) # dummy outside plot area
        end
        savefig(filename)
    end

    #=
    println("Plotting Θ0 + Ψ")
    plot(xlabel=L"x = \log a", ylabel=L"\Theta_0 + \Psi")
    ynum  = Θl.(co,x,k,0) .+ Ψ.(co,x,k)
    yanal = -maximum(abs.(ynum)) * cos.(k*c*η.(co,x)/√(3))
    plot!(x, ynum)
    plot!(x, yanal, color=:gray, linewidth=0.5)
    savefig("plots/Theta_plus_Psi.pdf")
    =#
end

end
