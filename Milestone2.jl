module Milestone2

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings
using Printf
using QuadGK

co_H_He_reion  = ΛCDM()
co_H_He_reioff = ΛCDM(z_reion_H=NaN)
co_H_reion     = ΛCDM(Yp=0)
co_H_reioff    = ΛCDM(Yp=0, z_reion_H=NaN)
co             = co_H_He_reion # default
x = range(-10, 0, length=8000)

# TODO: gather into one common function?
xswi = time_switch_Peebles(co)
xdec = time_decoupling(co) # TODO: compute from dg = 0
xrec = time_recombination(co)
xre1 = time_reionization_H(co)
xre2 = time_reionization_He(co)
shor = sound_horizon(co, xdec) / Gpc
println("Saha -> Peebles (Xe = 0.999):  ", format_time_variations(co, xswi))
println("Decoupling (max(g)):           ", format_time_variations(co, xdec))
println("Recombination (Xe = 0.1):      ", format_time_variations(co, xrec))
println("H  reionization (z_reion_H):   ", format_time_variations(co, xre1))
println("He reionization (z_reion_He):  ", format_time_variations(co, xre2))
println("Sound horizon at decoupling:   s = $shor Gpc, k = 2*π/s = $(2*π/shor) / Gpc")
println("Free electron fraction: today: ", Xe(co, 0))

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
    #scatter!( [time_switch_Peebles(co)],      [log10(0.999)], markersize=2, markerstrokewidth=0, color=:black, label=nothing)
    #annotate!([time_switch_Peebles(co)+0.15], [log10(0.999)+0.50], [(L"\textrm{Saha $\rightarrow$ Peebles}")])
    #annotate!([time_switch_Peebles(co)+0.10], [log10(0.999)+0.25], [(L"\textrm{($X_e = 0.999$)}")])

    # Annotate reionization on/off and recombination
    annotate!([-1, -1], [-0.2, -3.3], [("reionization"), ("reionizatioff")])
    annotate!([-8.5],   [-0.2], [("recombination(s)")])

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.5, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/free_electron_fraction_log.pdf")
end

if true || !isfile("plots/free_electron_fraction_linear.pdf")
    println("Plotting free electron fraction (linear)")

    plot(xlabel = L"x = \log a", ylabel = L"X_e", xlims=(x[1], x[end]), ylims=(-0.1, 1.3), yticks=-0.25:0.25:1.5, legend_position=:top, framestyle=:box)

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

    # Mark plateaus
    hline!([1+2*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{++}"])
    hline!([1+1*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}, \textrm{ He}^{+}"])
    hline!([1+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{H}^{+}"])
    hline!([0+0*co.Yp / (4*(1-co.Yp))], color = :gray, linestyle = :dash, z_order=:back, label=nothing) # label = [L"\textrm{fully ionized } \textrm{nothing}"])

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.5, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

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

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.5, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/optical_depth.pdf")
end

if true || !isfile("plots/visibility_function_linear.pdf")
    # TODO: inset reionization plot (simple example: http://www.breloff.com/images/juliacon/plotswithplots.slides.html#Inset/Floating-Subplots)
    println("Plotting visibility function (linear)")

    plot(xlabel = L"x = \log a", xlims=(x[1], x[end]), ylims=(-10, +10), legend_position=:topright)

    ys = [g.(co_H_reioff, x) / 1, dg.(co_H_reioff, x) / 10, d2g.(co_H_reioff, x) / 100,
          g.(co_H_He_reion,  x) / 1, dg.(co_H_He_reion,  x) / 10, d2g.(co_H_He_reion,  x) / 100]
    cs = [1  2  3  1  2  3]
    as = [0.3  0.3  0.3  1.0  1.0  1.0]
    Ls = [nothing  nothing  nothing  L"\tilde{g}\phantom{''}(x) \,/\, 1"  L"\tilde{g}'\phantom{'}(x) \,/\, 10"  L"\tilde{g}''(x) \,/\, 100"]
    zs = [:back  :back  :back  :back  :back  :back]

    plot!(x, ys, color=cs, alpha=as, z_order=zs, label=Ls)

    # Dummy plots to manually create legend
    hline!([-20], color=:black, linestyle=:solid, alpha=1.0, label=L"\textrm{H+He ($Y_p=0.24$), reionization}")
    hline!([-20], color=:black, linestyle=:solid, alpha=0.3, label=L"\textrm{H\phantom{+He} ($Y_p=0.00$), reionizatioff}")

    annotate!([-7.5], [9.0], [text(L"\int_{-20}^0 \tilde{g}(x) dx = %$(round(quadgk(x -> g(co_H_He_reion, x), -20, 0, rtol=1e-6)[1], digits=11))", color=:black)])
    annotate!([-7.5], [7.7], [text(L"\int_{-20}^0 \tilde{g}(x) dx = %$(round(quadgk(x -> g(co_H_reioff, x), -20, 0, rtol=1e-6)[1], digits=11))", color=:gray)])

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.25, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    # Zoom-in plot
    plot!(x, ys, color=cs, alpha=as, z_order=zs, label=nothing, xlims=(-3, -1), ylims=(-0.2, +0.2), subplot=2, inset = (1, bbox(0.09, 0.44, 0.3, 0.5, :right)))

    savefig("plots/visibility_function_linear.pdf")
end

if true || !isfile("plots/visibility_function_log.pdf")
    # TODO: inset reionization plot (simple example: http://www.breloff.com/images/juliacon/plotswithplots.slides.html#Inset/Floating-Subplots)
    println("Plotting visibility function (logarithmic)")

    plot(xlabel = L"x = \log a", ylabel = L"\log \tilde{g}", xlims=(x[1], x[end]), ylims=(-15, +5), legend_position=:topright)

    plot!(x, log10.(abs.(g.(co_H_He_reion, x))), color=1, alpha=1.0, label=nothing)
    plot!(x, log10.(abs.(g.(co_H_reioff, x))), color=1, alpha=0.3, label=nothing)

    # Mark event times
    vline!([xswi, xdec, xrec, xre1, xre2], linewidth=0.25, alpha=0.5, color=:black, linestyle=:dash, z_order=:back, label=nothing)

    savefig("plots/visibility_function_log.pdf")
end

end 
