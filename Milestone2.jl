module Milestone2

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings
using Printf

#co = ΛCDM(Ωb0=0.05, Ωc0=0.45, Neff=0, h=0.7)
co_H_He_reion  = ΛCDM()
co_H_He_reioff = ΛCDM(z_reion_H=NaN)
co_H_reion     = ΛCDM(Yp=0)
co_H_reioff    = ΛCDM(Yp=0, z_reion_H=NaN)
co             = co_H_He_reion # default
x = range(-10, 0, length=8000)

# TODO: gather into one common function?
xswi  = time_switch_Peebles(co)
xlss  = time_last_scattering_surface(co)
xrec  = time_recombination(co)
xdec  = (xlss + xrec) / 2
xre1  = time_reionization_H(co)
xre2  = time_reionization_He(co)
@printf("Reionization of Hydrogen:            x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xre1,  a(xre1),  z(xre1),  η(co, xre1)  / Gyr, t(co, xre1)  / Gyr)
@printf("Reionization of Helium:              x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xre2,  a(xre2),  z(xre2),  η(co, xre2)  / Gyr, t(co, xre2)  / Gyr)
@printf("Saha -> Peebles switch (Xe = 0.999): x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xswi,  a(xswi),  z(xswi),  η(co, xswi)  / Gyr, t(co, xswi)  / Gyr)
@printf("Last scattering surface (τ = 1):     x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xlss,  a(xlss),  z(xlss),  η(co, xlss)  / Gyr, t(co, xlss)  / Gyr)
@printf("Recombination (Xe = 0.1):            x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xrec,  a(xrec),  z(xrec),  η(co, xrec)  / Gyr, t(co, xrec)  / Gyr)
@printf("Decoupling (average(LSS, rec)):      x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xdec,  a(xdec),  z(xdec),  η(co, xdec)  / Gyr, t(co, xdec)  / Gyr)
println("Corresponding sound horizon:         $(sound_horizon(co, xdec) / Gpc) Gpc")
println("Freeze-out free electron fraction:   $(Xe(co, 0))")

if true || !isfile("plots/free_electron_fraction_log.pdf")
    println("Plotting free electron fraction (logarithmic)")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} X_e", xlims=(x[1], x[end]), ylims=(-4, 1.0), legendcolumns=2, legend_position=:topright)

    plot!(x, log10.(Xe.(co_H_reion,     x)),           linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization on
    plot!(x, log10.(Xe.(co_H_reioff,    x)),           linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization off
    plot!(x, log10.(Xe.(co_H_He_reion,  x)),           linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization on
    plot!(x, log10.(Xe.(co_H_He_reioff, x)),           linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization off
    plot!(x, log10.(Xe_Saha_H.(co_H_reioff, x)),       linestyle=:dash,  color=2, label=nothing) # Saha,         H,    reionization off
    plot!(x, log10.(Xe_Saha_H_He.(co_H_He_reioff, x)), linestyle=:dash,  color=1, label=nothing) # Saha,         H+He, reionization off # TODO: extend for x > -7

    # Dummy plots to manually create legend
    hline!([-10], color=:black, linestyle=:solid, label=L"\textrm{Saha+Peebles}")
    hline!([-10], color=1,      linestyle=:solid, label=L"\textrm{H+He ($Y_p=0.24$)}")
    hline!([-10], color=:black, linestyle=:dash,  label=L"\textrm{Saha}")
    hline!([-10], color=2,      linestyle=:solid, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$)}")

    # Mark and annotate Saha -> Peebles transition point
    scatter!( [time_switch_Peebles(co)],       [log10(0.999)], markersize=2, markerstrokewidth=0, color=:black, label=nothing)
    annotate!([time_switch_Peebles(co)+0.15], [log10(0.999)+0.50], [(L"\textrm{Saha $\rightarrow$ Peebles}")])
    annotate!([time_switch_Peebles(co)+0.10], [log10(0.999)+0.25], [(L"\textrm{($X_e = 0.999$)}")])

    # Annotate reionization on/off and recombination
    annotate!([-1, -1], [-0.2, -3.3], [("reionization"), ("reionizatioff")])
    annotate!([-8.5],   [-0.2], [("recombination(s)")])

    savefig("plots/free_electron_fraction_log.pdf")
end

if true || !isfile("plots/free_electron_fraction_linear.pdf")
    println("Plotting free electron fraction (linear)")

    plot(xlabel = L"x = \log a", ylabel = L"X_e", xlims=(x[1], x[end]), ylims=(-0.1, 1.3), yticks=-0.25:0.25:1.5, legend_position=:top, framestyle=:box)
    #plot!(x[x.<-7], Xe_Saha_H_He.(co, x[x.<-7]), label="Saha equation")
    #plot!(x, Xe.(co, x), label="Saha & Peebles equation")

    plot!(x, Xe.(co_H_reion,     x),           linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization on
    plot!(x, Xe.(co_H_reioff,    x),           linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization off
    plot!(x, Xe.(co_H_He_reion,  x),           linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization on
    plot!(x, Xe.(co_H_He_reioff, x),           linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization off
    plot!(x, Xe_Saha_H.(co_H_reioff, x),       linestyle=:dash,  color=2, label=nothing) # Saha,         H,    reionization off
    plot!(x, Xe_Saha_H_He.(co_H_He_reioff, x), linestyle=:dash,  color=1, label=nothing) # Saha,         H+He, reionization off # TODO: extend for x > -7

    # Annotate stages
    annotate!([-9.2], [1+2*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}^{++}$}}")])
    annotate!([-8.1], [1+1*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}^{+ }$}}")])
    annotate!([-7.4], [1+0*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}     $}}")])
    annotate!([-4.6], [0+0*co.Yp / (4*(1-co.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}  $, $\textrm{He}     $}}")])

    # Annotate reionization on/off
    annotate!([-1, -1], [0.95, 0.05], [("reionization"), ("reionizatioff")])
    annotate!([-8.5],   [0.95],       [("recombination(s)")])

    hline!([1+2*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{++}"])
    hline!([1+1*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{+}"])
    hline!([1+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}"])
    hline!([0+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{nothing}"])
    savefig("plots/free_electron_fraction_linear.pdf")
end

if true || !isfile("plots/optical_depth.pdf")
    println("Plotting optical depth")
    plot(xlabel = L"x = \log a", xlims=(x[1], x[end]), ylims=(-7.5, 3.5), legend_position=:topright)

    d2τpos(co, x) = d2τ(co, x) > 0 ? d2τ(co, x) : 1e-10 # skip values where d2τ < 0. increase resolution here?

    plot!(x, log10.(τ.(co_H_reioff, x)),      color=1, alpha=0.5, linestyle=:dash, label=nothing)
    plot!(x, log10.(-dτ.(co_H_reioff, x)),    color=2, alpha=0.5, linestyle=:dash, label=nothing)
    plot!(x, log10.(d2τpos.(co_H_reioff, x)), color=3, alpha=0.5, linestyle=:dash, label=nothing)

    plot!(x, log10.(τ.(co, x)),      color=1, linestyle=:solid, label=L"\log_{10} [+\tau\phantom{''}(x)]")
    plot!(x, log10.(-dτ.(co, x)),    color=2, linestyle=:solid, label=L"\log_{10} [-\tau'\phantom{'}(x)]")
    plot!(x, log10.(d2τpos.(co, x)), color=3, linestyle=:solid, label=L"\log_{10} [+\tau''(x) > 0]")

    # Dummy plots to manually create legend
    hline!([-10], color=:black, linestyle=:solid, alpha=1.0, label=L"\textrm{H+He ($Y_p=0.24$), reionization}")
    hline!([-10], color=:black, linestyle=:dash,  alpha=0.5, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$), reionizatioff}")

    savefig("plots/optical_depth.pdf")
end

if true || !isfile("plots/visibility_function.pdf")
    # TODO: inset reionization plot (simple example: http://www.breloff.com/images/juliacon/plotswithplots.slides.html#Inset/Floating-Subplots)
    println("Plotting visibility function")
    plot(xlabel = L"x = \log a", xlims=(x[1], x[end]), ylims=(-10, +10), legend_position=:topright)

    ys = [g.(co_H_He_reioff, x) / 1, dg.(co_H_He_reioff, x) / 10, d2g.(co_H_He_reioff, x) / 100,
          g.(co_H_He_reion,  x) / 1, dg.(co_H_He_reion,  x) / 10, d2g.(co_H_He_reion,  x) / 100]
    cs = [1  2  3  1  2  3]
    as = [0.3  0.3  0.3  1.0  1.0  1.0]
    Ls = [nothing  nothing  nothing  L"g\phantom{''}(x) \,/\, 1"  L"g'\phantom{'}(x) \,/\, 10"  L"g''(x) \,/\, 100"]
    zs = [:back  :back  :back  :back  :back  :back]

    plot!(x, ys, color=cs, alpha=as, z_order=zs, label=Ls)

    # Dummy plots to manually create legend
    hline!([-20], color=:black, linestyle=:solid, alpha=1.0, label=L"\textrm{H+He ($Y_p=0.24$), reionization}")
    hline!([-20], color=:black, linestyle=:solid, alpha=0.3, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$), reionizatioff}")

    #plot!(x, [], xlims=(-3, -1), ylims=(-0.2, +0.2), xticks=[-3,-2,-1], subplot=2, inset = (1, bbox(0.09, 0.5, 0.3, 0.5, :right)), label=nothing)
    plot!(x, ys, color=cs, alpha=as, z_order=zs, label=nothing, xlims=(-3, -1), ylims=(-0.2, +0.2), xticks=[-3,-2,-1], subplot=2, inset = (1, bbox(0.09, 0.43, 0.3, 0.5, :right)))

    savefig("plots/visibility_function.pdf")
end

end 
