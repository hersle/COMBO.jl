module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM(h=0.7, Neff=0, Ωb0=0.05, Ωc0=0.45, Yp=0, z_reion_H=NaN) # TODO: change back
k = 0.15 / Mpc
xtce = time_tight_coupling(co, k)
println("Tight coupling end (k = $(k*Mpc)/Mpc): ", format_time_variations(co, xtce))

if !isfile("plots/time_tight_coupling.pdf")
    println("Plotting tight coupling end time")
    log10ks = range(-3, +3, length=100)
    ks = 10 .^ log10ks / Mpc
    plot(xlabel=L"k / \textrm{Mpc}", ylabel=L"x")
    plot!(log10ks, time_tight_coupling.(co, ks))
    vline!([log10(0.15)])
    savefig("plots/time_tight_coupling.pdf")
end

if true || !isfile("plots/overdensity.pdf") || !isfile("plots/velocity.pdf") || !isfile("plots/potential.pdf") || !isfile("plots/temperature_fluctuation.pdf")
    k = 0.1 / Mpc
    x = range(-18, 0, length=5000)
    println("Using k = $k / Mpc")

    println("Plotting temperature fluctuation")
    plot(xlabel=L"x = \log a", ylabel=L"\Theta")
    # TODO: wrong sign with l=1?
    for l in 0:3
        plot!(x, Θl.(co, x, k, l), label=L"\Theta_%$(l)")
    end
    vline!([time_tight_coupling(co, k)], color=:gray, linestyle=:dash)
    savefig("plots/temperature_fluctuation.pdf")

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
