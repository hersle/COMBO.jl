module Milestone2

include("Cosmology.jl")

using .Cosmology
using Plots
using LaTeXStrings

co = ΛCDM()
x = range(-10, 0, length=10000)

if true || !isfile("plots/free_electron_fraction.pdf")
    println("Plotting free electron fraction")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} X_e", ylims=(-4, 0.5), legend_position=:topright)
    plot!(x, log10.(Xe.(co, x)), label="Saha and Peebles equation")
    plot!(x, log10.(Xe_Saha.(co, x)), linestyle=:dash, label="Saha equation")
    vline!([time_switch_Peebles(co)], linestyle=:dash, color=:gray, label = L"\textrm{Saha} \rightarrow \textrm{Peebles}")
    savefig("plots/free_electron_fraction.pdf")
end

if true || !isfile("plots/optical_depth.pdf")
    println("Plotting optical depth")
    plot(xlabel = L"x = \log a", legend_position=:topright)
    y = d2τ.(co, x[10:end-10])
    plot!(x, log10.(τ.(co, x)), label = L"\log_{10} [\tau(x)]")
    plot!(x, log10.(-dτ.(co, x)), label = L"\log_{10} [-\tau'(x)]")
    plot!(x, log10.(d2τ.(co, x)), label = L"\log_{10} [\tau''(x)]")
    savefig("plots/optical_depth.pdf")
end

if true || !isfile("plots/visibility_function.pdf")
    println("Plotting visibility function")
    plot(xlabel = L"x = \log a", ylabel = L"g", legend_position=:topright)
    plot!(x, g.(co, x), label = L"g(x)")
    plot!(x, dg.(co, x) / 10, label = L"g'(x) / 10") # TODO: handle endpoints
    plot!(x, d2g.(co, x) / 100, label = L"g''(x) / 10^2")
    savefig("plots/visibility_function.pdf")
end

end
