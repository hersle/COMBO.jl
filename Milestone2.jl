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
    y = d2τ.(co, x[10:end-10])
    println(y)
    println(minimum(y))
    println(maximum(y))
    plot!(x, log10.(τ.(co, x)))
    plot!(x[2:end-1], log10.(-dτ.(co, x[2:end-1])))
    plot!(x[10:end-10], log10.(d2τ.(co, x[10:end-10])))
    savefig("plots/optical_depth.pdf")
end

if true || !isfile("plots/visibility_function.pdf")
    println("Plotting visibility function")
    plot(xlabel = L"x = \log a", ylabel = L"g")
    plot!(x, g.(co, x))
    plot!(x[2:end-1], dg.(co, x[2:end-1])) # TODO: handle endpoints
    plot!(x[2:end-1], d2g.(co, x[2:end-1]))
    savefig("plots/visibility_function.pdf")
end

end
