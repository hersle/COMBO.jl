module Milestone1

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants

par = Parameters()
bg = Background(par)
x = range(-15, 5, length=400)
xrm = x_equality_rm(par)
xmΛ = x_equality_mΛ(par)
xacc = x_acceleration(par)
x1, x2, x3, x4 = minimum(x), xrm, xmΛ, maximum(x)
x0 = 0.0 # today

println("Ωr = Ωm:     ", format_time_variations(bg, xrm))
println("d²a/dt² = 0: ", format_time_variations(bg, xacc))
println("Ωm = ΩΛ:     ", format_time_variations(bg, xmΛ))
println("Today:       ", format_time_variations(bg, x0))

# Conformal Hubble parameter
if true
    println("Plotting plots/conformal_Hubble.pdf")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \mathcal{H} \,/\, (100\,\mathrm{km/s/Mpc}) \Big]", legend_position = :topleft, ylims = (-1, +7), yticks = -1:1:+7)
    plot!(x, @. log10(par.H0     / (100*km/Mpc) * √(par.Ωr0) * a(x)^(-1  )); linestyle = :dash,  label = "radiation-dominated")
    plot!(x, @. log10(par.H0     / (100*km/Mpc) * √(par.Ωm0) * a(x)^(-1/2)); linestyle = :dash,  label = "matter-dominated")
    plot!(x, @. log10(par.H0     / (100*km/Mpc) * √(par.ΩΛ0) * a(x)^(+1  )); linestyle = :dash,  label = "cosmological constant-dominated")
    plot!(x, @. log10(aH(par, x) / (100*km/Mpc)                          ); linestyle = :solid, label = "general case", color = :black)
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    vline!([xacc], z_order = :back, color = :gray, linestyle = :dot, label = nothing)
    annotate!([-2.5], [6.5], [(L"x_\textrm{acc} = %$(round(xacc, digits=2))", :gray)])
    savefig("plots/conformal_Hubble.pdf")
end

# conformal Hubble parameter 1st derivative
if true
    println("Plotting plots/conformal_Hubble_derivative1.pdf")
    plot(xlabel = L"x = \log a", ylabel = L"\frac{1}{\mathcal{H}} \frac{\mathrm{d}\mathcal{H}}{\mathrm{d} x}", legend_position = :topleft)
    plot!(x, x -> -1;                   linestyle = :dash,  label = "radiation-dominated")
    plot!(x, x -> -1/2;                 linestyle = :dash,  label = "matter-dominated")
    plot!(x, x -> +1;                   linestyle = :dash,  label = "cosmological constant-dominated")
    plot!(x, @. aH′(par, x) / aH(par, x); linestyle = :solid, label = "general case", color = :black)
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    vline!([xacc], z_order = :back, color = :gray, linestyle = :dot, label = nothing)
    savefig("plots/conformal_Hubble_derivative1.pdf")
end

# Conformal Hubble parameter 2nd derivative
if true
    println("Plotting plots/conformal_Hubble_derivative2.pdf")
    plot(xlabel = L"x = \log a", ylabel = L"\frac{1}{\mathcal{H}} \frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d} x^2}", legend_positions = :topleft, yticks = 0:0.25:1.5, ylims = (0, 1.5))
    plot!(x, x -> ( 1)^2;                linestyle = :dash,  label = "radiation-dominated")
    plot!(x, x -> (-1/2)^2;              linestyle = :dash,  label = "matter-dominated")
    plot!(x, x -> (-1)^2;                linestyle = :dash,  label = "cosmological constant-dominated", color = 1) # same color as radiation
    plot!(x, @. aH′′(par, x) / aH(par, x); linestyle = :solid, label = "general case", color = :black)
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    vline!([xacc], z_order = :back, color = :gray, linestyle = :dot, label = nothing)
    savefig("plots/conformal_Hubble_derivative2.pdf")
end

# Product of conformal time and conformal Hubble parameter
if true
    println("Plotting plots/eta_H.pdf")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \eta \mathcal{H} \Big]", legend_position = :topleft)
    plot!(x, @. log10(η(bg,x) * aH(par, x)); linestyle = :solid, color = :black, label = "general case")
    plot!(x, log10.(η_radiation_matter_dominated.(par, x) .* aH.(par, x)); linestyle = :dash, color = 1, label = "radiation-matter-dominated")
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    savefig("plots/eta_H.pdf")
end

# Cosmic and conformal time
if true
    println("Plotting plots/times.pdf")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \{t, \eta\} / \mathrm{Gyr} \Big]", legend_position = :bottomright)
    plot!(x, log10.(η.(bg, x)                             / Gyr); linestyle = :solid, color = 0, label = L"\eta \,\, \textrm{(general)}")
    plot!(x, log10.(η_radiation_matter_dominated.(par, x) / Gyr); linestyle = :dash,  color = 0, label = L"\eta \,\, \textrm{(radiation-matter-dominated)}")
    plot!(x, log10.(t.(bg, x)                             / Gyr); linestyle = :solid, color = 1, label = L"t \,\, \textrm{(general)}")
    plot!(x, log10.(t_radiation_matter_dominated.(par, x) / Gyr); linestyle = :dash,  color = 1, label = L"t \,\, \textrm{(radiation-matter-dominated}")
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    savefig("plots/times.pdf")
end

# Density parameters
if true
    println("Plotting plots/density_parameters.pdf")
    plot(xlabel = L"x = \log a", ylabel = L"\Omega_i", legend_position = (0.05, 0.6), ylims=(-0.05, +1.3))
    plot!(x, Ωr.(par, x); label = L"\Omega_r")
    plot!(x, Ωm.(par, x); label = L"\Omega_m")
    plot!(x, Ωk.(par, x); label = L"\Omega_k = %$(round(Int, Ωk(par, 0.0)))")
    plot!(x, ΩΛ.(par, x); label = L"\Omega_\Lambda")
    plot!(x, Ω.(par, x);  label = L"\sum_s \Omega_s = %$(round(Int, Ω(par, 0.0)))")
    plot!([xrm, xrm], [-0.05, 1.2]; z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    plot!([xmΛ, xmΛ], [-0.05, 1.2]; z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    annotate!([xrm], [1.25], [(L"x_\textrm{eq}^{rm} = %$(round(xrm; digits=2))", :gray)])
    annotate!([xmΛ], [1.25], [(L"x_\textrm{eq}^{m\Lambda} = %$(round(xmΛ; digits=2))", :gray)])
    annotate!([(x1+x2)/2, (x1+x2)/2], [1.14, 1.07], ["radiation", "domination"])
    annotate!([(x2+x3)/2, (x2+x3)/2], [1.14, 1.07], ["matter", "domination"])
    annotate!([(x3+x4)/2, (x3+x4)/2], [1.14, 1.07], ["Λ", "domination"])
    savefig("plots/density_parameters.pdf")
end

end
