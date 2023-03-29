module Milestone2

include("Cosmology.jl")
include("Constants.jl")

using .Cosmology
using .Constants
using Plots # TODO: common plotting settings. include("Plot.jl") with settings, or something?
using LaTeXStrings
using Printf

#co = ΛCDM(Ωb0=0.05, Ωc0=0.45, Neff=0, h=0.7)
co_H_He = ΛCDM()
co_H = ΛCDM(Yp=0)
#co = co_H_He # the default
co = ΛCDM()
x = range(-10, 0, length=10000)

# TODO: gather into one common function?
xswi  = time_switch_Peebles(co)
xlss  = time_last_scattering_surface(co)
xrec  = time_recombination(co)
xdec  = (xlss + xrec) / 2
xre1  = time_reionization_H(co)
xre2 = time_reionization_He(co)
@printf("Reionization of Hydrogen:           x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xre1,  a(xre1),  z(xre1),  η(co, xre1)  / Gyr, t(co, xre1)  / Gyr)
@printf("Reionization of Helium:             x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xre2,  a(xre2),  z(xre2),  η(co, xre2)  / Gyr, t(co, xre2)  / Gyr)
@printf("Saha -> Peebles switch (Xe = 0.99): x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xswi,  a(xswi),  z(xswi),  η(co, xswi)  / Gyr, t(co, xswi)  / Gyr)
@printf("Last scattering surface (τ = 1):    x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xlss,  a(xlss),  z(xlss),  η(co, xlss)  / Gyr, t(co, xlss)  / Gyr)
@printf("Recombination (Xe = 0.1):           x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xrec,  a(xrec),  z(xrec),  η(co, xrec)  / Gyr, t(co, xrec)  / Gyr)
@printf("Decoupling (average(LSS, rec)):     x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xdec,  a(xdec),  z(xdec),  η(co, xdec)  / Gyr, t(co, xdec)  / Gyr)
println("Corresponding sound horizon:        $(sound_horizon(co, xdec) / Gpc) Gpc")
println("Freeze-out free electron fraction:  $(Xe(co, 0; reionization=false))")

if !isfile("plots/free_electron_fraction_log.pdf")
    println("Plotting free electron fraction (logarithmic)")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} X_e", xlims=(x[1], x[end]), ylims=(-4, 1.0), legendcolumns=2, legend_position=:topright)

    plot!(x, log10.(Xe.(co_H, x; reionization=true)),    linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization on
    plot!(x, log10.(Xe.(co_H, x; reionization=false)),   linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization off
    plot!(x, log10.(Xe.(co, x; reionization=true)),      linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization on
    plot!(x, log10.(Xe.(co, x; reionization=false)),     linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization off
    plot!(x, log10.(Xe_Saha_H.(co_H, x)),                linestyle=:dash,  color=2, label=nothing) # Saha,         H,    reionization off
    plot!(x[x.<-7], log10.(Xe_Saha_H_He.(co, x[x.<-7])), linestyle=:dash,  color=1, label=nothing) # Saha,         H+He, reionization off # TODO: extend for x > -7

    # Dummy plots to manually create legend
    hline!([-10], color=:black, linestyle=:solid, label=L"\textrm{Saha+Peebles}")
    hline!([-10], color=1,      linestyle=:solid, label=L"\textrm{$Y_p=0.24$ (H+He)}")
    hline!([-10], color=:black, linestyle=:dash,  label=L"\textrm{Saha}")
    hline!([-10], color=2,      linestyle=:solid, label=L"\textrm{$Y_p=0.00$ (H)}")

    # Mark and annotate Saha -> Peebles transition point
    scatter!([time_switch_Peebles(co)], [log10(0.99)], markersize=2, markerstrokewidth=0, color=:black, label=nothing)
    annotate!([time_switch_Peebles(co)+0.15], [log10(0.99)+0.50], [(L"\textrm{Saha $\rightarrow$ Peebles}")])
    annotate!([time_switch_Peebles(co)+0.1], [log10(0.99)+0.25], [(L"\textrm{($X_e = 0.99$)}")])

    # Annotate reionization on/off and recombination
    annotate!([-1, -1], [-0.2, -3.3], [("reionization on"), ("reionization off")])
    annotate!([-8.5],   [-0.2], [("recombination(s)")])

    savefig("plots/free_electron_fraction_log.pdf")
end

if !isfile("plots/free_electron_fraction_linear.pdf")
    println("Plotting free electron fraction (linear)")

    plot(xlabel = L"x = \log a", ylabel = L"X_e", xlims=(x[1], x[end]), ylims=(-0.1, 1.3), yticks=-0.2:0.2:1.6, legend_position=:top, framestyle=:box)
    #plot!(x[x.<-7], Xe_Saha_H_He.(co, x[x.<-7]), label="Saha equation")
    #plot!(x, Xe.(co, x), label="Saha & Peebles equation")

    plot!(x,        Xe.(co_H, x; reionization=true),    linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization on
    plot!(x,        Xe.(co_H, x; reionization=false),   linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization off
    plot!(x,        Xe.(co, x; reionization=true),      linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization on
    plot!(x,        Xe.(co, x; reionization=false),     linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization off
    plot!(x,        Xe_Saha_H.(co_H, x),                linestyle=:dash,  color=2, label=nothing) # Saha,         H,    reionization off
    plot!(x[x.<-7], Xe_Saha_H_He.(co, x[x.<-7]),        linestyle=:dash,  color=1, label=nothing) # Saha,         H+He, reionization off # TODO: extend for x > -7

    # Annotate stages
    annotate!([-9.2], [1+2*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}^{++}$}}")])
    annotate!([-8.1], [1+1*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}^{+ }$}}")])
    annotate!([-7.4], [1+0*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}     $}}")])
    annotate!([-4.6], [0+0*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}  $, $\textrm{He}     $}}")])

    # Annotate reionization on/off
    annotate!([-1, -1], [0.95, 0.05], [("reionization on"), ("reionization off")])
    annotate!([-8.5],   [0.95],       [("recombination(s)")])

    hline!([1+2*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{++}"])
    hline!([1+1*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{+}"])
    hline!([1+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}"])
    hline!([0+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{nothing}"])
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
    # TODO: inset reionization plot (simple example: http://www.breloff.com/images/juliacon/plotswithplots.slides.html#Inset/Floating-Subplots)
    println("Plotting visibility function")
    plot(xlabel = L"x = \log a", ylabel = L"g", legend_position=:topright)
    plot!(x, g.(co, x), label = L"g(x)")
    plot!(x, dg.(co, x) / 10, label = L"g'(x) / 10") # TODO: handle endpoints
    plot!(x, d2g.(co, x) / 100, label = L"g''(x) / 10^2")
    savefig("plots/visibility_function.pdf")
end

end 
