module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM(h=0.7, Neff=0, Ωb0=0.05, Ωc0=0.45, Yp=0, z_reion_H=NaN) # TODO: change back
#co = ΛCDM()
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

if true || !isfile("plots/overdensity.pdf") || !isfile("plots/velocity.pdf") || !isfile("plots/potential.pdf") || !isfile("plots/temperature_fluctuation.pdf") # TODO: add more temperature fluctuations...
    settings = [
        ("plots/velocity.pdf",  L"\log_{10}|v|", [(Cosmology.i_vc, log10, :solid, L"v=v_c"), (Cosmology.i_vb, log10∘abs, :dash, L"v=v_b")]),
        ("plots/velocity2.pdf", L"\log_{10}|v|", [(Cosmology.i_Θl(1), Θ1 -> log10(abs(-3*Θ1)), :solid, L"v=v_\gamma=-3\Theta_1"), (Cosmology.i_Nl(1), N1 -> log10(abs(-3*N1)), :dash, L"v=v_\nu=-3\mathcal{N}_1")]),
        ("plots/overdensity.pdf",  L"\log_{10}|\delta|", [(Cosmology.i_δc, log10, :solid, L"\delta=\delta_c"), (Cosmology.i_δb, log10∘abs, :dash, L"\delta=\delta_b")]),
        ("plots/overdensity2.pdf", L"\log_{10}|\delta|", [(Cosmology.i_Θl(0), Θ0 -> log10(abs(+4*Θ0)), :solid, L"\delta=\delta_\gamma= 4 \Theta_0"), (Cosmology.i_Nl(0), N0 -> log10(abs(+4*N0)), :dash, L"\delta=\delta_\nu= 4\mathcal{N}_0")]),
        ("plots/Nl0.pdf", L"\mathcal{N}_0", [(Cosmology.i_Nl(0), identity, :solid, "")]),
        ("plots/Nl1.pdf", L"\mathcal{N}_1", [(Cosmology.i_Nl(1), identity, :solid, "")]),
        ("plots/Nl2.pdf", L"\mathcal{N}_2", [(Cosmology.i_Nl(2), identity, :solid, "")]),
        ("plots/ThetaPl0.pdf", L"\Theta^P_0", [(Cosmology.i_ΘPl(0), identity, :solid, "")]),
        ("plots/ThetaPl1.pdf", L"\Theta^P_1", [(Cosmology.i_ΘPl(1), identity, :solid, "")]),
        ("plots/ThetaPl2.pdf", L"\Theta^P_2", [(Cosmology.i_ΘPl(2), identity, :solid, "")]),
        ("plots/Thetal0.pdf", L"\Theta_0", [(Cosmology.i_Θl(0), identity, :solid, "")]),
        ("plots/Thetal1.pdf", L"\Theta_1", [(Cosmology.i_Θl(1), identity, :solid, "")]),
        ("plots/Thetal2.pdf", L"\Theta_2", [(Cosmology.i_Θl(2), identity, :solid, "")]),
        ("plots/potential.pdf", L"\Phi", [(Cosmology.i_Φ, identity, :solid, "")])
    ]

    # pre-compute callable splines once and for all (index and call as y1s[i_k][i_qty](x))
    y1s = [[x -> spline(x, k) for spline in Cosmology.perturbations_splines(co; tight=false)] for k in ks] # 2D (x, k) splines for full system
    y2s = [Cosmology.perturbations_mode(co, k; tight=true ) for k in ks] # 1D (x) splines for each "exact" k for tight+full system
    y3s = [Cosmology.perturbations_mode(co, k; tight=false) for k in ks] # 1D (x) splines for each "exact" k for       full system

    for (filename, ylabel, iqty_func_linestyle_labels) in settings
        println("Plotting $filename")
        plot(xlabel=L"x = \log a", ylabel=ylabel, xlims=(-20, 0), xticks=-25:5:5, legend_position=:topleft)
        for (i, k) in enumerate(ks)
            x = extendx(Cosmology.splinex(y3s[i][1]), 7) # plot with extra points between every spline point for more smoothness
            label = L"k = %$(k*Mpc) / \textrm{Mpc}"
            for (i_qty, func, linestyle, _) in iqty_func_linestyle_labels
                plot!(x, func.(y1s[i][i_qty].(x)), alpha=1/8, linewidth=3.0, color=i, linestyle=linestyle, label=nothing)
                plot!(x, func.(y2s[i][i_qty].(x)), alpha=2/8, linewidth=2.0, color=i, linestyle=linestyle, label=nothing)
                plot!(x, func.(y3s[i][i_qty].(x)), alpha=8/8, linewidth=1.0, color=i, linestyle=linestyle, label=label)
                label = nothing # only label each k-value once
            end
            vline!([time_tight_coupling(co, k)], color=:gray, linestyle=:dash; linewidth=0.5, label=nothing)
        end
        # add quantity (if more than one so ambiguous) to legend with a black dummy plot
        if length(iqty_func_linestyle_labels) > 1
            for (i_qty, func, linestyle, label) in iqty_func_linestyle_labels
                vline!([-21], linewidth=1.0, linestyle=linestyle, color=:black, label=label) # dummy outside plot area
            end
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
