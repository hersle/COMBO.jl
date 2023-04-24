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

function plot_quantity_with_all_methods!(co, x, k, i_qty::Integer; label=nothing, kwargs...)
    plot!(x, Cosmology.perturbations_splines(co)[i_qty].(x, k),             alpha=2/6, linewidth=1.5, label=nothing; kwargs...)
    plot!(x, Cosmology.perturbations_mode(co, k, 6; tight=true )[i_qty](x), alpha=4/6, linewidth=1.0, label=nothing; kwargs...)
    plot!(x, Cosmology.perturbations_mode(co, k, 6; tight=false)[i_qty](x), alpha=6/6, linewidth=0.5, label=label; kwargs...)
end

if true || !isfile("plots/overdensity.pdf") || !isfile("plots/velocity.pdf") || !isfile("plots/potential.pdf") || !isfile("plots/temperature_fluctuation.pdf") # TODO: add more temperature fluctuations...
    println("Plotting temperature fluctuations")
    x = range(-20.0, 0.0, length=5000)
    for l in 0:3
        plot(xlabel=L"x = \log a", ylabel=L"\Theta_%$(l)")
        for (i, k) in enumerate(ks)
            plot_quantity_with_all_methods!(co, x, k, Cosmology.i_Θl(0); label=L"k = %$(k*Mpc) / \textrm{Mpc}", color=i)
            vline!([time_tight_coupling(co, k)], color=i, linestyle=:dash; label=nothing)
        end
        savefig("plots/temperature_fluctuation_l$l.pdf")
    end

    println("Plotting Θ0 + Ψ")
    plot(xlabel=L"x = \log a", ylabel=L"\Theta_0 + \Psi")
    ynum  = Θl.(co,x,k,0) .+ Ψ.(co,x,k)
    yanal = -maximum(abs.(ynum)) * cos.(k*c*η.(co,x)/√(3))
    plot!(x, ynum)
    plot!(x, yanal, color=:gray, linewidth=0.5)
    savefig("plots/Theta_plus_Psi.pdf")

    # return

    println("Plotting overdensity")
    plot(xlabel=L"x = \log a", ylabel=L"\log_{10} \delta")
    plot!(x, log10.(δc.(co, x, k)))
    plot!(x, log10.(abs.(δb.(co, x, k))))
    savefig("plots/overdensity.pdf")

    println("Plotting velocity")
    plot(xlabel=L"x = \log a", ylabel=L"\log_{10} v")
    plot!(x, log10.(vc.(co, x, k)))
    plot!(x, log10.(abs.(vb.(co, x, k))))
    savefig("plots/velocity.pdf")

    println("Plotting potential")
    plot(xlabel=L"x = \log a", ylabel=L"\Phi")
    plot!(x, Φ.(co, x, k))
    plot!(x, Φ.(co, x, k))
    savefig("plots/potential.pdf")
end

end
