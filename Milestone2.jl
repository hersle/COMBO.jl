module Milestone2

include("Cosmology.jl")
include("Constants.jl")

using .Cosmology
using .Constants
using Plots # TODO: common plotting settings. include("Plot.jl") with settings, or something?
using LaTeXStrings
using Printf

#co = ΛCDM(Ωb0=0.05, Ωc0=0.45, Neff=0, h=0.7)
co = ΛCDM()
x = range(-10, 0, length=10000)

# TODO: gather into one common function?
xswi = time_switch_Peebles(co)
xlss = time_last_scattering_surface(co)
xrec = time_recombination(co)
xdec = (xlss + xrec) / 2
@printf("Saha -> Peebles switch (Xe ) 0.99): x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xswi,  a(xswi),  z(xswi),  η(co, xswi)  / Gyr, t(co, xswi)  / Gyr)
@printf("Last scattering surface (τ = 1):    x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xlss,  a(xlss),  z(xlss),  η(co, xlss)  / Gyr, t(co, xlss)  / Gyr)
@printf("Recombination (Xe = 0.1):           x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xrec,  a(xrec),  z(xrec),  η(co, xrec)  / Gyr, t(co, xrec)  / Gyr)
@printf("Decoupling (average(LSS, rec)):     x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xdec,  a(xdec),  z(xdec),  η(co, xdec)  / Gyr, t(co, xdec)  / Gyr)
println("Corresponding sound horizon:        $(sound_horizon(co, xdec) / Gpc) Gpc")
println("Freeze-out free electron fraction:  $(Xe(co, 0; reionization=false))")

if true || !isfile("plots/free_electron_fraction_log.pdf") || !isfile("plots/free_electron_fraction_linear.pdf")
    println("Plotting free electron fraction")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} X_e", ylims=(-4, 0.5), legend_position=:topright)
    plot!(x, log10.(Xe.(co, x)), label="Saha & Peebles equation")
    plot!(x, log10.(Xe_Saha_H.(co, x)), linestyle=:dash, label="Saha equation")
    vline!([time_switch_Peebles(co)], linestyle=:dash, color=:gray, label = L"\textrm{Saha} \rightarrow \textrm{Peebles}")
    savefig("plots/free_electron_fraction_log.pdf")

    plot(xlabel = L"x = \log a", ylabel = L"X_e", ylims=(-0.1, 1.5), legend_position=:top)
    plot!(x[x.<-7], Xe_Saha_H_He.(co, x[x.<-7]), label="Saha equation")
    plot!(x, Xe.(co, x), label="Saha & Peebles equation")
    hline!([1+2*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{++}"])
    hline!([1+1*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{+}"])
    hline!([1+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, label = [L"\textrm{fully ionized } \textrm{H}^{+}"])
    hline!([0+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, label = [L"\textrm{fully ionized } \textrm{nothing}"])
    savefig("plots/free_electron_fraction_linear.pdf")
end

if true || !isfile("plots/optical_depth.pdf")
    println("Plotting optical depth")
    plot(xlabel = L"x = \log a", legend_position=:topright)
    plot!(x, log10.(τ.(co, x)), label = L"\log_{10} [\tau(x)]")
    plot!(x, log10.(-dτ.(co, x)), label = L"\log_{10} [-\tau'(x)]")
    y(x) = d2τ(co, x) > 0 ? log10(d2τ(co, x)) : NaN # skip values where d2τ < 0. increase resolution here?
    plot!(x, y.(x), label = L"\log_{10} [\tau''(x) > 0]") # TODO: is the little kink wrong?
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
