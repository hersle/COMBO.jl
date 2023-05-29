module Milestone2

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings
using Printf
using QuadGK

rec_H_He_reion  = Recombination(Background(Parameters())) # hydrogen, helium and reionization (default)
rec_H_He_reioff = Recombination(Background(Parameters(z_reion_H=NaN))) # no reionization
rec_H_reion     = Recombination(Background(Parameters(Yp=0))) # no helium
rec_H_reioff    = Recombination(Background(Parameters(Yp=0, z_reion_H=NaN))) # no reionization and helium
rec_Saha        = Recombination(Background(Parameters()), xswitch=0.0) # never switch to Peebles
rec             = rec_H_He_reion # default

xswi = x_switch_Peebles(rec.bg.par)
xdec = x_decoupling(rec)
xrec = x_recombination(rec)
xre1 = x_reionization_H(rec.bg.par)
xre2 = x_reionization_He(rec.bg.par)
xdec_Saha = x_decoupling(rec_Saha)
xrec_Saha = x_recombination(rec_Saha)
shor = s(rec, xdec) / Gpc
println("Saha -> Peebles (Xe = 0.999):         ", format_time_variations(rec.bg, xswi))
println("Decoupling (max(g)):                  ", format_time_variations(rec.bg, xdec))
println("Decoupling (max(g)) (Saha only):      ", format_time_variations(rec_Saha.bg, xdec_Saha))
println("Recombination (Xe = 0.1):             ", format_time_variations(rec.bg, xrec))
println("Recombination (Xe = 0.1) (Saha only): ", format_time_variations(rec_Saha.bg, xrec_Saha))
println("H  reionization (z_reion_H):          ", format_time_variations(rec.bg, xre1))
println("He reionization (z_reion_He):         ", format_time_variations(rec.bg, xre2))
println("Sound horizon at decoupling:          s = $shor Gpc, k = 2*π/s = $(2*π/shor) / Gpc")
println("Free electron fraction: today:        ", Xe(rec, 0))

if true || !isfile("plots/free_electron_fraction_log.pdf")
    println("Plotting free electron fraction (logarithmic)")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} X_e", xlims=(-10, 0), ylims=(-4, 1.0), legendcolumns=2, legend_position=:topright)

    x = range(-20, 0, length=1000)
    plot!(x, log10.(Xe.(rec_H_reion,     x)),                  linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization on
    plot!(x, log10.(Xe.(rec_H_reioff,    x)),                  linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization off
    plot!(x, log10.(Xe.(rec_H_He_reion,  x)),                  linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization on
    plot!(x, log10.(Xe.(rec_H_He_reioff, x)),                  linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization off
    plot!(x, log10.(Xe_Saha_H.(rec_H_reioff.bg.par, x)),       linestyle=:dash,  color=2, label=nothing) # Saha,         H,    reionization off
    plot!(x, log10.(Xe_Saha_H_He.(rec_H_He_reioff.bg.par, x)), linestyle=:dash,  color=1, label=nothing) # Saha,         H+He, reionization off # TODO: extend for x > -7

    # Dummy plots to manually create legend
    hline!([-10], color=:black, linestyle=:solid, label=L"\textrm{Saha+Peebles}")
    hline!([-10], color=1,      linestyle=:solid, label=L"\textrm{H+He ($Y_p=0.24$)}")
    hline!([-10], color=:black, linestyle=:dash,  label=L"\textrm{Saha}")
    hline!([-10], color=2,      linestyle=:solid, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$)}")

    # Annotate reionization on/off and recombination
    annotate!([-1, -1], [-0.2, -3.3], [("reionization"), ("reionizatioff")])
    annotate!([-8.5],   [-0.2], [("recombination(s)")])

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.5, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/free_electron_fraction_log.pdf")
end

if true || !isfile("plots/free_electron_fraction_linear.pdf")
    println("Plotting free electron fraction (linear)")

    plot(xlabel = L"x = \log a", ylabel = L"X_e", xlims=(-10, 0), ylims=(-0.1, 1.3), yticks=-0.25:0.25:1.5, legend_position=:top, framestyle=:box)

    x = range(-20, 0, length=1000)
    plot!(x, Xe.(rec_H_reion,     x),                  linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization on
    plot!(x, Xe.(rec_H_reioff,    x),                  linestyle=:solid, color=2, label=nothing) # Saha+Peebles, H,    reionization off
    plot!(x, Xe.(rec_H_He_reion,  x),                  linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization on
    plot!(x, Xe.(rec_H_He_reioff, x),                  linestyle=:solid, color=1, label=nothing) # Saha+Peebles, H+He, reionization off
    plot!(x, Xe_Saha_H.(rec_H_reioff.bg.par, x),       linestyle=:dash,  color=2, label=nothing) # Saha,         H,    reionization off
    plot!(x, Xe_Saha_H_He.(rec_H_He_reioff.bg.par, x), linestyle=:dash,  color=1, label=nothing) # Saha,         H+He, reionization off # TODO: extend for x > -7

    # Annotate stages
    annotate!([-9.2], [1+2*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}^{++}$}}")])
    annotate!([-8.1], [1+1*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}^{+ }$}}")])
    annotate!([-7.4], [1+0*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}^+$, $\textrm{He}     $}}")])
    annotate!([-4.6], [0+0*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))+0.03], [(L"\small{\textrm{$\textrm{H}  $, $\textrm{He}     $}}")])

    # Annotate reionization on/off
    annotate!([-1, -1], [0.95, 0.05], [("reionization"), ("reionizatioff")])
    annotate!([-8.5],   [0.95],       [("recombination(s)")])

    # Mark plateaus
    hline!([1+2*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{++}"])
    hline!([1+1*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{+}"])
    hline!([1+0*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}"])
    hline!([0+0*rec.bg.par.Yp / (4*(1-rec.bg.par.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{nothing}"])

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.5, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/free_electron_fraction_linear.pdf")
end

if true || !isfile("plots/optical_depth.pdf")
    println("Plotting optical depth")

    plot(xlabel = L"x = \log a", xlims=(-10, 0), ylims=(-7.5, 3.5), legend_position=:topright)

    x = Cosmology.integration_points(rec.τ)
    τ′′pos(rec, x) = τ′′(rec, x) > 0 ? τ′′(rec, x) : 1e-10 # "positive" τ′′: skip values where τ′′ < 0. increase resolution here?

    plot!(x, log10.(+τ.(rec_H_reioff, x)),      color=1, alpha=0.5, linestyle=:dash, label=nothing)
    plot!(x, log10.(-τ′.(rec_H_reioff, x)),     color=2, alpha=0.5, linestyle=:dash, label=nothing)
    plot!(x, log10.(+τ′′pos.(rec_H_reioff, x)), color=3, alpha=0.5, linestyle=:dash, label=nothing)

    plot!(x, log10.(+τ.(rec, x)),      color=1, linestyle=:solid, label=L"\log_{10} [+\tau\phantom{''}(x)]")
    plot!(x, log10.(-τ′.(rec, x)),     color=2, linestyle=:solid, label=L"\log_{10} [-\tau'\phantom{'}(x)]")
    plot!(x, log10.(+τ′′pos.(rec, x)), color=3, linestyle=:solid, label=L"\log_{10} [+\tau''(x) > 0]")

    τreion, Δτreion = 0.054, 0.007
    #hline!([log10(τreion)], ribbon=([log10(τreion)-log10(τreion-Δτreion)], [log10(τreion+Δτreion)-log10(τreion)]), color=:gray, label=L"\tau_\textrm{reion} (Planck 2018)")
    hline!([log10(τreion)], color=:gray, alpha=0.5, label=L"\log_{10} [\tau_\textrm{reion}] \textrm{ (Planck 2018)}")

    # Dummy plots to manually create legend
    hline!([-10], color=:black, linestyle=:solid, alpha=1.0, label=L"\textrm{H+He ($Y_p=0.24$), reionization}")
    hline!([-10], color=:black, linestyle=:dash,  alpha=0.5, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$), reionizatioff}")

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.5, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/optical_depth.pdf")
end

if true || !isfile("plots/visibility_function_linear.pdf")
    println("Plotting visibility function (linear)")

    plot(xlabel = L"x = \log a", xlims=(-10, 0), ylims=(-10, +10), legend_position=:topright)

    x = Cosmology.integration_points(rec.τ)
    ys = [g.(rec_H_reioff, x) / 1, g′.(rec_H_reioff, x) / 10, g′′.(rec_H_reioff, x) / 100,
          g.(rec_H_He_reion,  x) / 1, g′.(rec_H_He_reion,  x) / 10, g′′.(rec_H_He_reion,  x) / 100]
    cs = [1  2  3  1  2  3]
    as = [0.3  0.3  0.3  1.0  1.0  1.0]
    Ls = [nothing  nothing  nothing  L"\tilde{g}\phantom{''}(x) \,/\, 1"  L"\tilde{g}'\phantom{'}(x) \,/\, 10"  L"\tilde{g}''(x) \,/\, 100"]
    zs = [:back  :back  :back  :back  :back  :back]

    plot!(x, ys, color=cs, alpha=as, z_order=zs, label=Ls)

    # Dummy plots to manually create legend
    hline!([-20], color=:black, linestyle=:solid, alpha=1.0, label=L"\textrm{H+He ($Y_p=0.24$), reionization}")
    hline!([-20], color=:black, linestyle=:solid, alpha=0.3, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$), reionizatioff}")

    annotate!([-7.5], [9.0], [text(L"\int_{-20}^0 \tilde{g}(x) dx = %$(round(quadgk(x -> g(rec_H_He_reion, x), -20, 0, rtol=1e-6)[1], digits=11))", color=:black)])
    annotate!([-7.5], [7.7], [text(L"\int_{-20}^0 \tilde{g}(x) dx = %$(round(quadgk(x -> g(rec_H_reioff, x),   -20, 0, rtol=1e-6)[1], digits=11))", color=:gray)])

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.25, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    # Zoom-in reionization plot
    # (simple example: http://www.breloff.com/images/juliacon/plotswithplots.slides.html#Inset/Floating-Subplots)
    x = range(-3, -1, length=400) # re-plot with many points to resolve narrow peaks
    ys = [g.(rec_H_reioff, x) / 1, g′.(rec_H_reioff, x) / 10, g′′.(rec_H_reioff, x) / 100,
          g.(rec_H_He_reion,  x) / 1, g′.(rec_H_He_reion,  x) / 10, g′′.(rec_H_He_reion,  x) / 100]
    plot!(x, ys, color=cs, alpha=as, z_order=zs, label=nothing, xlims=extrema(x), ylims=(-0.2, +0.2), subplot=2, inset = (1, bbox(0.09, 0.44, 0.3, 0.5, :right)))

    savefig("plots/visibility_function_linear.pdf")
end

if true || !isfile("plots/visibility_function_log.pdf")
    println("Plotting visibility function (logarithmic)")

    plot(xlabel = L"x = \log a", ylabel = L"\log \tilde{g}", xlims=(-10, 0), ylims=(-15, +5), legend_position=:topright)

    x = Cosmology.integration_points(rec.τ)
    plot!(x, log10.(abs.(g.(rec_H_He_reion, x))), color=1, alpha=1.0, label=nothing)
    plot!(x, log10.(abs.(g.(rec_H_reioff, x))), color=1, alpha=0.3, label=nothing)

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.25, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/visibility_function_log.pdf")
end

end 
