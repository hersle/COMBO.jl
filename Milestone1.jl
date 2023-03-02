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
using Distributions
using Printf
Plots.__init__() # workaround with sysimage: https://github.com/JuliaLang/PackageCompiler.jl/issues/786
pgfplotsx()
default(
    minorticks = 10,
    labelfontsize = 12, # default: 11
    legendfontsize = 11, # default: 8
    annotationfontsize = 11, # default: 14
    tickfontsize = 10, # default: 8
    #titlefontsize = 14, # default: 14
    legend_font_halign = :left,
)

# a hack to "rasterize" a scatter plot
# TODO: plotting scatter points as rects lowers file size?
function scatterheatmaps!(xss, yss, colors, labels, xlims, ylims; nbins=50, kwargs...)
    xgrid = range(xlims...; length=nbins+1)
    ygrid = range(ylims...; length=nbins+1)
    dx = xgrid[2] - xgrid[1]
    dy = ygrid[2] - ygrid[1]
    zgrid = zeros(Int, nbins, nbins)
    for (i, xs, ys) in zip(1:length(xss), xss, yss)
        for (x, y) in zip(xs, ys)
            ix = min(floor(Int, (x - xgrid[1]) / dx) + 1, nbins)
            iy = min(floor(Int, (y - ygrid[1]) / dy) + 1, nbins)
            zgrid[iy,ix] = i
        end
    end
    xgrid = (xgrid[1:end-1] .+ xgrid[2:end]) ./ 2
    ygrid = (ygrid[1:end-1] .+ ygrid[2:end]) ./ 2
    #zgrid = zgrid .!= 0
    colorgrad = cgrad([:white, colors...], alpha=1.0, categorical=false)
    #histogram2d!(xs, ys; bins=(range(xlims...; length=nbins), range(ylims...; length=nbins)), color = colorgrad, colorbar = :none, kwargs...)
    heatmap!(xgrid, ygrid, zgrid; alpha = 1.0, color = colorgrad, colorbar = :none, kwargs...)
    for (color, label) in zip(colors, labels)
        scatter!([xgrid[1]-1], [ygrid[1]-1]; color = color, markerstrokewidth=0, shape = :rect, label = label) # dummy plot with label
    end
    # TODO: hack a label into the legend with a ghost scatter plot?
end

co = ΛCDM()
x = range(-15, 5, length=400)
xrm = r_m_equality(co)
xmΛ = m_Λ_equality(co)
xacc = acceleration_onset(co)
x1, x2, x3, x4 = minimum(x), xrm, xmΛ, maximum(x)
x0 = 0.0

@printf("Ωr = Ωm:     x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xrm,  a(xrm),  z(xrm),  η(co, xrm)  / Gyr, t(co, xrm)  / Gyr)
@printf("d²a/dt² = 0: x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xacc, a(xacc), z(xacc), η(co, xacc) / Gyr, t(co, xacc) / Gyr)
@printf("Ωm = ΩΛ:     x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", xmΛ,  a(xmΛ),  z(xmΛ),  η(co, xmΛ)  / Gyr, t(co, xmΛ)  / Gyr)
@printf("Today:       x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr\n", x0,   a(x0),   z(x0),   η(co, x0)   / Gyr, t(co, x0)   / Gyr)


if !isdir("plots")
    mkdir("plots")
end

# Conformal Hubble parameter
if true || !isfile("plots/conformal_Hubble.pdf")
    println("Plotting conformal Hubble parameter")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10}\Big[ \mathcal{H} \,/\, (100\,\mathrm{km/s/Mpc})\Big]", legend_position = :topleft, ylims = (-1, +7), yticks = -1:1:+7)
    plot!(x, @. log10(co.H0     / (100*km/Mpc) * √(co.Ωr0) * a(x)^(-1  )); linestyle = :dash,  label = "radiation-dominated")
    plot!(x, @. log10(co.H0     / (100*km/Mpc) * √(co.Ωm0) * a(x)^(-1/2)); linestyle = :dash,  label = "matter-dominated")
    plot!(x, @. log10(co.H0     / (100*km/Mpc) * √(co.ΩΛ0) * a(x)^(+1  )); linestyle = :dash,  label = "cosmological constant-dominated")
    plot!(x, @. log10(aH(co, x) / (100*km/Mpc)                          ); linestyle = :solid, label = "general case", color = :black)
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    savefig("plots/conformal_Hubble.pdf")
end

# conformal Hubble parameter 1st derivative
if true || !isfile("plots/conformal_Hubble_derivative1.pdf")
    println("Plotting conformal Hubble 1st derivative")
    plot(xlabel = L"x = \log a", ylabel = L"\frac{1}{\mathcal{H}} \frac{\mathrm{d}\mathcal{H}}{\mathrm{d} x}", legend_position = :topleft)
    plot!(x, x -> -1;                   linestyle = :dash,  label = "radiation-dominated")
    plot!(x, x -> -1/2;                 linestyle = :dash,  label = "matter-dominated")
    plot!(x, x -> +1;                   linestyle = :dash,  label = "cosmological constant-dominated")
    plot!(x, @. daH(co, x) / aH(co, x); linestyle = :solid, label = "general case", color = :black)
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    savefig("plots/conformal_Hubble_derivative1.pdf")
end

# Conformal Hubble parameter 2nd derivative
if true || !isfile("plots/conformal_Hubble_derivative2.pdf")
    println("Plotting conformal Hubble 2nd derivative")
    plot(xlabel = L"x = \log a", ylabel = L"\frac{1}{\mathcal{H}} \frac{\mathrm{d}^2\mathcal{H}}{\mathrm{d} x^2}", legend_positions = :topleft, yticks = 0:0.25:1.5, ylims = (0, 1.5))
    plot!(x, x -> ( 1)^2;                linestyle = :dash,  label = "radiation-dominated")
    plot!(x, x -> (-1/2)^2;              linestyle = :dash,  label = "matter-dominated")
    plot!(x, x -> (-1)^2;                linestyle = :dash,  label = "cosmological constant-dominated", color = 1) # same color as radiation
    plot!(x, @. d2aH(co, x) / aH(co, x); linestyle = :solid, label = "general case", color = :black)
    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    savefig("plots/conformal_Hubble_derivative2.pdf")
end

# Product of conformal time and conformal Hubble parameter
if true || !isfile("plots/eta_H.pdf")
    println("Plotting η * aH")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \eta \mathcal{H} \Big]", legend_position = :topleft)

    #plot!(x, x -> log10(1);                      linestyle = :dash,  label = "radiation-dominated")
    plot!(x, @. log10(η(co, x) * aH(co, x)); linestyle = :solid, color = :black, label = "general case")

    aeq_anal = co.Ωr0 / co.Ωm0
    η_anal = @. 2 / (co.H0 * √(co.Ωm0)) * (√(a(x) + aeq_anal) - √(aeq_anal))
    aH_anal = @. a(x) * co.H0 * √(co.Ωr0/a(x)^4 + co.Ωm0/a(x)^3)
    η_aH_anal = @. η_anal * aH_anal
    plot!(x, log10.(η_aH_anal); linestyle = :dash, color = 1, label = "radiation-matter universe")

    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)

    savefig("plots/eta_H.pdf")
end

# Cosmic and conformal time
if true || !isfile("plots/times.pdf")
    println("Plotting cosmic and conformal times")
    plot(xlabel = L"x = \log a", ylabel = L"\log_{10} \Big[ \{t, \eta\} / \mathrm{Gyr} \Big]", legend_position = :bottomright)

    # in a radiation-matter-only universe
    aeq_anal = co.Ωr0 / co.Ωm0
    η_anal = @. 2 / (    co.H0 * √(co.Ωm0)) * (√(a(x) + aeq_anal) - √(aeq_anal))
    t_anal = @. 2 / (3 * co.H0 * √(co.Ωm0)) * (√(a(x) + aeq_anal) * (a(x) - 2*aeq_anal) + 2*aeq_anal^(3/2))

    plot!(x, log10.(η.(co, x) / Gyr); linestyle = :solid, color = 0, label = L"\eta \,\, \textrm{(general)}")
    plot!(x, log10.(η_anal    / Gyr); linestyle = :dash,  color = 0, label = L"\eta \,\, \textrm{(radiation-matter universe)}")

    plot!(x, log10.(t.(co, x) / Gyr); linestyle = :solid, color = 1, label = L"t \,\, \textrm{(general)}")
    plot!(x, log10.(t_anal    / Gyr); linestyle = :dash,  color = 1, label = L"t \,\, \textrm{(radiation-matter universe)}")

    vline!([xrm, xmΛ], z_order = :back, color = :gray, linestyle = :dash, label = nothing)

    savefig("plots/times.pdf")
end

# Density parameters
if true || !isfile("plots/density_parameters.pdf")
    println("Plotting density parameters")
    plot(xlabel = L"x = \log a", ylabel = L"\Omega_i", legend_position = (0.05, 0.6), ylims=(-0.05, +1.3))
    plot!(x, Ωr.(co, x); label = L"\Omega_r")
    plot!(x, Ωm.(co, x); label = L"\Omega_m")
    plot!(x, Ωk.(co, x); label = L"\Omega_k = %$(round(Int, Ωk(co, 0.0)))")
    plot!(x, ΩΛ.(co, x); label = L"\Omega_\Lambda")
    plot!(x, Ω.(co, x);  label = L"\sum_s \Omega_s = %$(round(Int, Ω(co, 0.0)))")
    plot!([xrm, xrm], [-0.05, 1.2]; z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    plot!([xmΛ, xmΛ], [-0.05, 1.2]; z_order = :back, color = :gray, linestyle = :dash, label = nothing)
    annotate!([xrm], [1.25], [(L"x_{r=m} = %$(round(xrm; digits=2))", :gray)])
    annotate!([xmΛ], [1.25], [(L"x_{m=\Lambda} = %$(round(xmΛ; digits=2))", :gray)])
    annotate!([(x1+x2)/2, (x1+x2)/2], [1.14, 1.07], ["radiation", "domination"])
    annotate!([(x2+x3)/2, (x2+x3)/2], [1.14, 1.07], ["matter", "domination"])
    annotate!([(x3+x4)/2, (x3+x4)/2], [1.14, 1.07], ["Λ", "domination"])
    savefig("plots/density_parameters.pdf")
end

# Supernova data
if true || !isfile("plots/supernova_distance.pdf")
    #co = ΛCDM(h=0.7, Ωc0=0.30-0.05, Ωk0=0.00)
    plot(xlabel = L"\log_{10} \Big[ 1+z \Big]", ylabel = L"\log_{10} \Big[ d_L \,/\, \mathrm{Gpc} \Big]", xlims=(0, 0.4), ylims=(-1.5, 1.0), legend_position = :topleft)

    x2 = range(-1, 0, length=400)
    plot!(log10.(Cosmology.z.(x2).+1), log10.(Cosmology.dL.(co, x2) ./ Gpc); label = "prediction")

    data = readdlm("data/supernovadata.txt", comments=true)
    zobs, dL, σdL = data[:,1], data[:,2], data[:,3]
    err_lo = log10.(dL) - log10.(dL-σdL)
    err_hi = log10.(dL+σdL) - log10.(dL)
    scatter!(log10.(zobs.+1), log10.(dL); markercolor = :black, yerror = (err_lo, err_hi), markersize=2, label = "supernova observations")

    savefig("plots/supernova_distance.pdf")
end

# Supernova MCMC fits
# Inspiration: "A theoretician's analysis of the supernova data ..." (https://arxiv.org/abs/astro-ph/0212573)
if true || !isfile("plots/supernova_omegas.pdf") || !isfile("plots/supernova_hubble.pdf")
    println("Plotting Ωm0, ΩΛ from MCMC analysis of supernova data")

    data = readdlm("data/supernovadata.txt", comments=true)
    N_obs, _ = size(data)
    z_obs, dL_obs, σdL_obs = data[:,1], data[:,2], data[:,3]
    x_obs = @. -log(z_obs + 1)

    function logLfunc(params::Vector{Float64})
        h, Ωm0, Ωk0 = params[1], params[2], params[3]
        co = ΛCDM(h=h, Ωb0=0.05, Ωc0=Ωm0-0.05, Ωk0=Ωk0, Neff=0)
        if Cosmology.is_fucked(co)
            return -Inf # so set L = 0 (or log(L) = -∞, or χ2 = ∞)
        else
            dL_mod = Cosmology.dL.(co, x_obs) / Gpc
            return -1/2 * sum((dL_mod .- dL_obs).^2 ./ σdL_obs.^2) # L = exp(-χ2/2)
        end
    end

    hbounds = (0.5, 1.5)
    Ωm0bounds = (0.0, 1.0)
    Ωk0bounds = (-1.0, +1.0)
    nparams = 3 # h, Ωm0, Ωk0
    nchains = 5
    nsamples = 10000 # per chain
    params, logL = MetropolisHastings(logLfunc, [hbounds, Ωm0bounds, Ωk0bounds], nsamples, nchains; burnin=1000)
    h, Ωm0, Ωk0, χ2 = params[:,1], params[:,2], params[:,3], -2 * logL

    # compute corresponding ΩΛ values by reconstructing the cosmologies (this is computationally cheap)
    ΩΛ0 = [ΛCDM(h=h[i], Ωb0=0.05, Ωc0=Ωm0[i]-0.05, Ωk0=Ωk0[i], Neff=0).ΩΛ0 for i in 1:length(h)]
    ΩΛ0bounds = (0.0, 1.2) # just for later plotting

    # best fit
    best_index = argmax(logL)
    best_h, best_Ωm0, best_Ωk0, best_ΩΛ0, best_χ2 = h[best_index], Ωm0[best_index], Ωk0[best_index], ΩΛ0[best_index], χ2[best_index]
    println("Best fit (χ²/N = $(round(best_χ2/N_obs, digits=1))): h = $best_h, Ωm0 = $best_Ωm0, Ωk0 = $best_Ωk0")

    # compute 68% ("1σ") and 95% ("2σ") confidence regions (https://en.wikipedia.org/wiki/Confidence_region)
    # (we have 3D Gaussian, but only for a 1D Gaussian does this correspond to the probability of finding values within 1σ/2σ!)
    # (see e.g. https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    #  and      https://stats.stackexchange.com/questions/61273/68-confidence-level-in-multinormal-distributions
    #  and      https://math.stackexchange.com/questions/1357071/confidence-ellipse-for-a-2d-gaussian)
    confidence2 = cdf(Chisq(1), 2^2) # ≈ 95% (2σ away from mean of 1D Gaussian)
    confidence1 = cdf(Chisq(1), 1^2) # ≈ 68% (1σ away from mean of 1D Gaussian)
    filter2 = χ2 .- best_χ2 .< quantile(Chisq(nparams), confidence2) # ≈ 95% ("2σ" in a 1D Gaussian) confidence region
    filter1 = χ2 .- best_χ2 .< quantile(Chisq(nparams), confidence1) # ≈ 68% ("1σ" in a 1D Gaussian) confidence region
    h2, Ωm02, Ωk02, ΩΛ02 = h[filter2], Ωm0[filter2], Ωk0[filter2], ΩΛ0[filter2]
    h1, Ωm01, Ωk01, ΩΛ01 = h[filter1], Ωm0[filter1], Ωk0[filter1], ΩΛ0[filter1]

    plot(xlabel = L"\Omega_{m0}", ylabel = L"\Omega_{\Lambda}", size = (600, 600), xlims = ΩΛ0bounds, ylims = ΩΛ0bounds, xticks = range(ΩΛ0bounds..., step=0.1), yticks = range(ΩΛ0bounds..., step=0.1), legend_position = :topright)
    #scatter!(Ωm02, ΩΛ02; color = 1, markershape = :rect, markerstrokecolor = 1, markerstrokewidth = 0, markersize = 2.0, clip_mode = "individual", label = L"%$(round(confidence2*100; digits=1)) % \textrm{ confidence region}") # clip mode workaround to get line above scatter points: https://discourse.julialang.org/t/plots-jl-with-pgfplotsx-adds-series-in-the-wrong-order/85896
    #scatter!(Ωm01, ΩΛ01; color = 3, markershape = :rect, markerstrokecolor = 3, markerstrokewidth = 0, markersize = 2.0, clip_mode = "individual", label = L"%$(round(confidence1*100; digits=1)) % \textrm{ confidence region}")
    scatterheatmaps!([Ωm02, Ωm01], [ΩΛ02, ΩΛ01], [palette(:default)[1], :darkblue], [L"%$(round(confidence2*100; digits=1)) % \textrm{ confidence region}", L"%$(round(confidence1*100; digits=1)) % \textrm{ confidence region}"], ΩΛ0bounds, ΩΛ0bounds; nbins=120)


    # plot ΩΛ(Ωm0) for a few flat universes (should give ΩΛ ≈ 1 - Ωm0)
    Ωm0_flat = range(Ωm0bounds..., length=20)
    ΩΛ0_flat = [ΛCDM(h=best_h, Ωb0=0.05, Ωc0=Ωm0-0.05, Ωk0=0, Neff=0).ΩΛ0 for Ωm0 in Ωm0_flat]
    plot!(Ωm0_flat, ΩΛ0_flat; color = :black, marker = :circle, markersize = 2, label = "flat universes")

    scatter!([best_Ωm0], [best_ΩΛ0]; color = :red, markerstrokecolor = :red, markershape = :cross, markersize = 10, label = "our best fit")
    scatter!([0.317], [0.683]; color = :green, markerstrokecolor = :green, markershape = :cross, markersize = 10, label = "Planck 2018's best fit")
    #vline([best_Ωm0]; linestyle = :dash, color = :red, label = L"\textrm{our best fit}")
    #hline([best_ΩΛ0]; linestyle = :dash, color = :red)
    savefig("plots/supernova_omegas.pdf")

    # TODO: draw error ellipses?
    # see e.g. https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

    println("Plotting h from MCMC analysis of supernova data")

    nfit = fit(Normal, h) # fit Hubble parameters to normal distribution
    plot(xlabel = L"h = H_0 \,/\, 100\,\frac{\mathrm{km}}{\mathrm{s}\,\mathrm{Mpc}}", ylabel = L"P(h)", xlims = (0.65, 0.75), ylims = (0, 1.1 / √(2*π*var(nfit))), xticks = 0.65:0.01:0.75, yticks=0:10:100, legend_position = :topright)
    histogram!(h; normalize = true, linewidth = 0, label = L"%$(nchains) \times %$(nsamples) \textrm{ samples}")
    vline!([best_h]; linestyle = :dash, color = :red, label = L"\textrm{our best fit}")
    vline!([0.674]; linestyle = :dash, label = L"\textrm{Planck 2018's best fit}")
    plot!(h -> pdf(nfit, h); color = :black, label = L"N(\mu = %$(round(mean(nfit), digits=2)), \sigma = %$(round(std(nfit), digits=2)))")
    savefig("plots/supernova_hubble.pdf")
end

end
