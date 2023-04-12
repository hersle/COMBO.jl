module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM(Neff=0)
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
    k = 0.01 / Mpc
    x = range(-18, time_tight_coupling(co,k), length=200)
    println("Using k = $k / Mpc")

    println("Plotting overdensity")
    plot(xlabel=L"x = \log a", ylabel=L"\log_{10} \delta")
    plot!(x, log10.(δc.(co, x, k)))
    plot!(x, log10.(δb.(co, x, k)))
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

    println("Plotting temperature fluctuation")
    plot(xlabel=L"x = \log a", ylabel=L"\Theta")
    plot!(x, Θl.(co, x, k, 0))
    plot!(x, Θl.(co, x, k, 1)) # TODO: wrong sign?
    savefig("plots/temperature_fluctuation.pdf")
end

end
