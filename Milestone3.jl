module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM(h=0.7, Neff=0, Ωb0=0.05, Ωc0=0.45, Yp=0, z_reion_H=NaN) # TODO: change back
ks = [1e-3, 1e-2, 1e-1] / Mpc
#println("Tight coupling end (k = $(k*Mpc)/Mpc): ", format_time_variations(co, xtce))

function plot_quantity_with_all_methods!(co, x, k, i_qty::Integer; func=y->y, label=nothing, kwargs...)
    plot!(x, func.(Cosmology.perturbations_splines(co)[i_qty].(x, k)),             alpha=2/6, linewidth=1.5, label=nothing; kwargs...)
    plot!(x, func.(Cosmology.perturbations_mode(co, k, 6; tight=true )[i_qty](x)), alpha=4/6, linewidth=1.0, label=nothing; kwargs...)
    plot!(x, func.(Cosmology.perturbations_mode(co, k, 6; tight=false)[i_qty](x)), alpha=6/6, linewidth=0.5, label=label; kwargs...)
end

if true || !isfile("plots/overdensity.pdf") || !isfile("plots/velocity.pdf") || !isfile("plots/potential.pdf") || !isfile("plots/temperature_fluctuation.pdf") # TODO: add more temperature fluctuations...
    x = range(-20.0, 0.0, length=5000)

    settings = [
        ("plots/temperature_fluctuation_l0.pdf", L"\Theta_0", [(Cosmology.i_Θl(0), identity, :solid)]),
        ("plots/temperature_fluctuation_l1.pdf", L"\Theta_1", [(Cosmology.i_Θl(1), identity, :solid)]),
        ("plots/temperature_fluctuation_l2.pdf", L"\Theta_2", [(Cosmology.i_Θl(2), identity, :solid)]),
        ("plots/overdensity.pdf", L"\delta", [(Cosmology.i_δc, log10, :solid), (Cosmology.i_δb, log10∘abs, :dash)]),
        ("plots/velocity.pdf", L"v", [(Cosmology.i_vc, log10, :solid), (Cosmology.i_vb, log10∘abs, :dash)]),
        ("plots/potential.pdf", L"\Phi", [(Cosmology.i_Φ, identity, :solid)])
    ]

    for (filename, ylabel, iqty_func_linestyles) in settings
        println("Plotting ", filename)
        plot(xlabel=L"x = \log a", ylabel=ylabel, xticks=-25:5:5)
        for (i, k) in enumerate(ks)
            for (i_qty, func, linestyle) in iqty_func_linestyles
                plot_quantity_with_all_methods!(co, x, k, i_qty; func=func, label=L"k = %$(k*Mpc) / \textrm{Mpc}", color=i, linestyle=linestyle)
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
