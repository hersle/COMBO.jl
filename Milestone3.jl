module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM()
k = 0.15 / Mpc
xtce = time_tight_coupling(co, k)
println("Tight coupling end: ", format_time_variations(co, xtce))

if !isfile("plots/time_tight_coupling.pdf")
    println("Plotting tight coupling end time")
    log10ks = range(-3, +3, length=100)
    ks = 10 .^ log10ks / Mpc
    plot(xlabel=L"k / \textrm{Mpc}", ylabel=L"x")
    plot!(log10ks, time_tight_coupling.(co, ks))
    vline!([log10(0.15)])
    savefig("plots/time_tight_coupling.pdf")
end

#perturbations(co, k, 1e-18, 10)

end
