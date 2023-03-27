module Milestone2

include("Cosmology.jl")

using .Cosmology
using Plots
using LaTeXStrings

co = ΛCDM()
x = range(-15, 0, length=1000)

if true || !isfile("plots/free_electron_fraction.pdf")
    println("Plotting free electron fraction")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} X_e", ylims=(-4, 0.5))
    plot!(x, log10.(Xe.(co, x)), label="Saha and Peebles")
    plot!(x, log10.(Xe_Saha.(co, x)), linestyle=:dash, label="Saha")
    savefig("plots/free_electron_fraction.pdf")
end

if true || !isfile("plots/optical_depth.pdf")
    println("Plotting optical depth")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \tau")
    plot!(x, log10.(τ.(co, x)))
    savefig("plots/optical_depth.pdf")
end

end
