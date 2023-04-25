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
        ("plots/velocity.pdf", L"v", [(Cosmology.i_vc, log10, :solid), (Cosmology.i_vb, log10∘abs, :dash)]),
        ("plots/N0.pdf", L"N_0", [(Cosmology.i_Nl(0), identity, :solid)]),
        ("plots/N1.pdf", L"N_1", [(Cosmology.i_Nl(1), identity, :solid)]),
        ("plots/N2.pdf", L"N_2", [(Cosmology.i_Nl(2), identity, :solid)]),
        ("plots/ThetaPl0.pdf", L"\Theta_0", [(Cosmology.i_ΘPl(0), identity, :solid)]),
        ("plots/ThetaPl1.pdf", L"\Theta_1", [(Cosmology.i_ΘPl(1), identity, :solid)]),
        ("plots/ThetaPl2.pdf", L"\Theta_2", [(Cosmology.i_ΘPl(2), identity, :solid)]),
        ("plots/Thetal0.pdf", L"\Theta_0", [(Cosmology.i_Θl(0), identity, :solid)]),
        ("plots/Thetal1.pdf", L"\Theta_1", [(Cosmology.i_Θl(1), identity, :solid)]),
        ("plots/Thetal2.pdf", L"\Theta_2", [(Cosmology.i_Θl(2), identity, :solid)]),
        ("plots/overdensity.pdf", L"\delta", [(Cosmology.i_δc, log10, :solid), (Cosmology.i_δb, log10∘abs, :dash)]),
        ("plots/potential.pdf", L"\Phi", [(Cosmology.i_Φ, identity, :solid)])
    ]

    # pre-compute callable splines once and for all (index and call as y1s[i_k][i_qty](x))
    y1s = [[x -> spline(x, k) for spline in Cosmology.perturbations_splines(co; tight=false)] for k in ks] # 2D (x, k) splines for full system
    y2s = [Cosmology.perturbations_mode(co, k; tight=true ) for k in ks] # 1D (x) splines for each "exact" k for tight+full system
    y3s = [Cosmology.perturbations_mode(co, k; tight=false) for k in ks] # 1D (x) splines for each "exact" k for       full system

    for (filename, ylabel, iqty_func_linestyles) in settings
        println("Plotting ", filename)
        plot(xlabel=L"x = \log a", ylabel=ylabel, xticks=-25:5:5)
        for (i, k) in enumerate(ks)
            label = L"k = %$(k*Mpc) / \textrm{Mpc}"
            x = extendx(Cosmology.splinex(y3s[i][1]), 3) # plot with extra points between every spline point for more smoothness
            for (i_qty, func, linestyle) in iqty_func_linestyles
                plot!(x, func.(y1s[i][i_qty].(x)), alpha=2/6, linewidth=1.5, label=nothing, color=i, linestyle=linestyle)
                plot!(x, func.(y2s[i][i_qty].(x)), alpha=4/6, linewidth=1.0, label=nothing, color=i, linestyle=linestyle)
                plot!(x, func.(y3s[i][i_qty].(x)), alpha=6/6, linewidth=0.5, label=label,   color=i, linestyle=linestyle)
            end
            vline!([time_tight_coupling(co, k)], color=:gray, linestyle=:dash; linewidth=0.5, label=nothing)
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
