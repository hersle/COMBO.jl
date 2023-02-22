module Milestone1

include("Algorithms.jl")
include("Cosmology.jl")
include("Constants.jl")

using .Cosmology
using .Algorithms
using .Constants
using Plots
using LaTeXStrings
using DelimitedFiles
Plots.__init__() # workaround with sysimage: https://github.com/JuliaLang/PackageCompiler.jl/issues/786
pgfplotsx()
default(legend_font_halign = :left)

co = ΛCDM()
x = range(-15, 5, length=400)

if !isdir("plots")
    mkdir("plots")
end

# Conformal Hubble parameter
if !isfile("plots/conformal_Hubble.pdf")
    println("Plotting conformal Hubble parameter")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10}\Big[ \mathcal{H} \,/\, (100\,\mathrm{km/s/Mpc})\Big]", legend_position = :topright)
    plot!(x, @. log10(co.H0     / (100*km/Mpc) * √(co.Ωr0) * a(x)^(-1  )); linestyle = :dash,  label = "radiation-dominated")
    plot!(x, @. log10(co.H0     / (100*km/Mpc) * √(co.Ωm0) * a(x)^(-1/2)); linestyle = :dash,  label = "matter-dominated")
    plot!(x, @. log10(co.H0     / (100*km/Mpc) * √(co.ΩΛ)  * a(x)^(+1  )); linestyle = :dash,  label = "cosmological constant-dominated")
    plot!(x, @. log10(aH(co, x) / (100*km/Mpc)                          ); linestyle = :solid, label = "general case", color = :black)
    savefig("plots/conformal_Hubble.pdf")
end

# conformal Hubble parameter 1st derivative
if !isfile("plots/conformal_Hubble_derivative1.pdf")
    println("Plotting conformal Hubble 1st derivative")
    plot(xlabel = L"x = \log a", ylabel = L"\frac{1}{\mathcal{H}} \frac{\mathrm{d}\mathcal{H}}{\mathrm{d} x}", legend_position = :topleft)
    plot!(x, x -> -1;                   linestyle = :dash,  label = "radiation-dominated")
    plot!(x, x -> -1/2;                 linestyle = :dash,  label = "matter-dominated")
    plot!(x, x -> +1;                   linestyle = :dash,  label = "cosmological constant-dominated")
    plot!(x, @. daH(co, x) / aH(co, x); linestyle = :solid, label = "general case", color = :black)
    savefig("plots/conformal_Hubble_derivative1.pdf")
end

# Conformal Hubble parameter 2nd derivative
if !isfile("plots/conformal_Hubble_derivative2.pdf")
    println("Plotting conformal Hubble 2nd derivative")
    plot(xlabel = L"x = \log a", ylabel = L"\frac{1}{\mathcal{H}} \frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d} x^2}", legend_positions = :topleft)
    plot!(x, x -> ( 1)^2;                linestyle = :dash,  label = "radiation and cosmological constant-dominated")
    plot!(x, x -> (-1/2)^2;              linestyle = :dash,  label = "matter-dominated")
    plot!(x, @. d2aH(co, x) / aH(co, x); linestyle = :solid, label = "general case", color = :black)
    savefig("plots/conformal_Hubble_derivative2.pdf")
end

# Product of conformal time and conformal Hubble parameter
if !isfile("plots/eta_H.pdf")
    println("Plotting η * aH / c")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \eta \mathcal{H} / c \Big]", legend_position = :topleft)

    #plot!(x, x -> log10(1);                      linestyle = :dash,  label = "radiation-dominated")
    plot!(x, @. log10(η(co, x) * aH(co, x) / c); linestyle = :solid, label = "general case")

    aeq_anal = co.Ωr0 / co.Ωm0
    η_anal = @. 2*c / (co.H0 * √(co.Ωm0)) * (√(a(x) + aeq_anal) - √(aeq_anal))
    aH_anal = @. a(x) * co.H0 * √(co.Ωr0/a(x)^4 + co.Ωm0/a(x)^3)
    η_aH_c_anal = @. η_anal * aH_anal / c
    plot!(x, log10.(η_aH_c_anal); linestyle = :dash, label = "radiation- and matter domination")

    savefig("plots/eta_H.pdf")
end

# Cosmic and conformal time
if !isfile("plots/times.pdf")
    println("Plotting cosmic and conformal times")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \{t, \eta/c\} / \mathrm{Gyr} \Big]", legend_position = :bottomright)

    # in a radiation-matter-only universe
    aeq_anal = co.Ωr0 / co.Ωm0
    η_anal = @. 2*c / (co.H0 * √(co.Ωm0)) * (√(a(x) + aeq_anal) - √(aeq_anal))
    t_anal = @.   2 / (3 * co.H0 * √(co.Ωm0)) * (√(a(x) + aeq_anal) * (a(x) - 2*aeq_anal) + 2*aeq_anal^(3/2))

    plot!(x, log10.(η.(co, x) / c / Gyr); linestyle = :solid, color = 0, label = L"\eta \,\, \mathrm{(general)}")
    plot!(x, log10.(η_anal    / c / Gyr); linestyle = :dash,  color = 0, label = L"\eta \,\, \mathrm{(radiation-matter universe)}")

    plot!(x, log10.(t.(co, x) / Gyr); linestyle = :solid, color = 1, label = L"t \,\, \mathrm{(general)}")
    plot!(x, log10.(t_anal    / Gyr); linestyle = :dash,  color = 1, label = L"t \,\, \mathrm{(radiation-matter universe)}")
    savefig("plots/times.pdf")
end

# Density parameters
if !isfile("plots/density_parameters.pdf")
    println("Plotting density parameters")
    plot(xlabel = L"x = \log a", ylabel = L"\Omega_i", legend_position = (0.05, 0.6), ylims=(-0.05, +1.3))
    plot!(x, Ωr.(co, x); label = L"\Omega_r")
    plot!(x, Ωm.(co, x); label = L"\Omega_m")
    plot!(x, Ωk.(co, x); label = L"\Omega_k")
    plot!(x, ΩΛ.(co, x); label = L"\Omega_\Lambda")
    plot!(x, Ω.(co, x);  label = L"\sum_s \Omega_s")
    xrm = r_m_equality(co)
    xmΛ = m_Λ_equality(co)
    x1, x2, x3, x4 = minimum(x), xrm, xmΛ, maximum(x)
    plot!([xrm, xrm], [-0.05, 1.2]; z_order = :back, color = :gray, linestyle = :dash, label = L"\Omega_r = \Omega_m")
    plot!([xmΛ, xmΛ], [-0.05, 1.2]; z_order = :back, color = :gray, linestyle = :dash, label = L"\Omega_m = \Omega_\Lambda")
    annotate!([xrm, xmΛ], [1.25, 1.25], [L"x = %$(round(xrm; digits=1))", L"x = %$(round(xmΛ; digits=1))"])
    annotate!([(x1+x2)/2, (x1+x2)/2], [1.15, 1.05], ["radiation", "domination"])
    annotate!([(x2+x3)/2, (x2+x3)/2], [1.15, 1.05], ["matter", "domination"])
    annotate!([(x3+x4)/2, (x3+x4)/2], [1.15, 1.05], ["Λ", "domination"])
    savefig("plots/density_parameters.pdf")
end

# Supernova data
if !isfile("plots/supernova.pdf")
    #co = ΛCDM(h=0.7, Ωc0=0.31-0.05, Ωk0=-0.03)
    plot(xlabel = L"\log_{10} \Big[ 1+z \Big]", ylabel = L"\log_{10} \Big[ d_L \,/\, \mathrm{Gpc} \Big]", xlims=(0, 0.4), ylims=(-1.5, 1.5), legend_position = :topleft)

    x2 = range(-1, 0, length=400)
    plot!(log10.(Cosmology.z.(x2).+1), log10.(Cosmology.dL.(co, x2) ./ Gpc); label = "prediction")

    data = readdlm("data/supernovadata.txt", comments=true)
    z, dL, σdL = data[:,1], data[:,2], data[:,3]
    err_lo = log10.(dL) - log10.(dL-σdL)
    err_hi = log10.(dL+σdL) - log10.(dL)
    scatter!(log10.(z.+1), log10.(dL); yerror = (err_lo, err_hi), markersize=2, label = "supernova observations")

    savefig("plots/supernova.pdf")
end

# Supernova MCMC fits
if true || !isfile("plots/supernova_mcmc.pdf") || !isfile("plots/supernova_hubble.pdf")
    data = readdlm("data/supernovadata.txt", comments=true)
    N_obs, _ = size(data)
    z_obs, dL_obs, σdL_obs = data[:,1], data[:,2], data[:,3]
    x_obs = @. log(1 / (z_obs+1))

    function logL(params::Vector{Float64})
        h, Ωm0, Ωk0 = params # unpack
        co = ΛCDM(h=h, Ωb0=0.05, Ωc0=Ωm0-0.05, Ωk0=Ωk0)
        if isnan(Cosmology.dL.(co, maximum(x_obs))) # model does not extend far enough so that it can fit the data
            return -Inf # so set L = 0 (or log(L) = -∞, or χ2 = ∞)
        else
            dL_mod = Cosmology.dL.(co, x_obs) / Gpc
            return -1/2 * sum(@. (dL_mod - dL_obs)^2 / σdL_obs^2) # L = exp(-χ2/2) # TODO: optimize
        end
    end

    params, logLs = MetropolisHastings(logL, ([0.5, 0.0, -1.0], [1.5, 1.0, +1.0]), 2000; steps=[0.007, 0.05, 0.05])
    params = params[1:end-1000, :] # remove burn-in
    logLs  =  logLs[1:end-1000]
    h, Ωm0, Ωk0 = params[:,1], params[:,2], params[:,3]
    best_index = argmax(logLs)
    best_χ2 = -2 * logLs[best_index]
    best_h, best_Ωm0, best_Ωk0 = params[best_index,:]
    println("Best fit (χ²/N = $(round(best_χ2/N_obs, digits=1))): h = $best_h, Ωm0 = $best_Ωm0, Ωk0 = $best_Ωk0")

    # compute corresponding ΩΛ values by reconstructing the cosmologies
    ΩΛ = [ΛCDM(h=h[i], Ωb0=0.05, Ωc0=Ωm0[i]-0.05, Ωk0=Ωk0[i]).ΩΛ for i in 1:length(h)]
    plot(xlabel = L"\Omega_{m0}", ylabel = L"\Omega_{\Lambda}", xlims = (0, 1), ylims = (-1, +1))
    scatter!(Ωm0, ΩΛ; label = nothing)
    savefig("plots/supernova_mcmc.pdf")

    plot(xlabel = L"h = H_0 \,/\, 100\,\frac{\mathrm{km}}{\mathrm{s}\,\mathrm{Mpc}}", ylabel = L"P(h)",xlims = (0.65, 0.75))
    histogram!(h; normalize = true, label = nothing)
    savefig("plots/supernova_hubble.pdf")
end

end
