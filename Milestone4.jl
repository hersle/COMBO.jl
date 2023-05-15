module Milestone4

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings
using DelimitedFiles

par = Parameters()
bg = Background(par)
rec = Recombination(bg)
xrm = equality_rm(par)
krm = 1 / (c*η(bg,xrm))

function plot_Dl(series, polarization, filename)
    plot(xlabel=L"\log_{10} l", ylabel=L"C_l^\mathrm{%$(polarization)} \, T_{\gamma 0} \cdot \frac{l(l+1)}{2\pi} \,/\, (\mathrm{\mu} K)^2", legend_position=:topleft)

    for (i, (l, Dl, ΔDl, settings)) in enumerate(series)
        println("i=$i, l = $l")
        plot!(log10.(l), Dl; yerror=ΔDl, color=i, markerstrokecolor=i, markersize=1, marker=:circle, settings...)
    end

    savefig(filename)
end

function lgrid(; n1=10, n2=20, n3=150)
    return unique(Int.(round.(10 .^ vcat(
        range(0.0, 1.0, length=n1),
        range(1.0, 2.0, length=n2),
        range(2.0, 3.4, length=n3)
    ))))
end

function plot_Dl_against_Planck(type)
    @assert type in (:TT, :TE, :EE)

    l = lgrid()

    plot_Dl([
        (l, Dl(rec,l,type), nothing, Dict(:label => "Our ΛCDM prediction")),
        (read_Planck_data("data/Planck_Dl_$type.txt")..., Dict(:label => "Planck's measurements", :color => :black, :marker => :square, :markerstrokecolor => :black, :markersize => 0.75, :markerstrokewidth => 0.50, :seriestype => :scatter)),
    ], "$type", "plots/power_spectrum_CMB_$type.pdf")
end

function plot_Dl_varying_parameter(paramkey, values; labelfunc=val -> nothing, type=:TT)
    l = lgrid()
    series = []
    for (i, value) in enumerate(values)
        par = Parameters(; Dict(paramkey => value)...)
        bg = Background(par)
        rec = Recombination(bg)

        # * With 3 values, the idea is always to plot
        #   the Planck values[2] in "neutral" black,
        #   the lower values[1] in "colder" blue,
        #   the higher values[3] in "warmer" red
        # * With 2 values, the idea is to plot
        #   a   feature  enabled with values[2] in "neutral" black,
        #   the feature disabled with values[1] in "cold" blue
        color = [:blue, :black, :red][i]

        push!(series, (l, Dl(rec,l,type), nothing, Dict(:color => color, :marker => :none, :label => labelfunc(value))))
    end
    plot_Dl(series, "$type", "plots/power_spectrum_CMB_varying_$paramkey.pdf")
end

function read_Planck_data(filename)
    data = readdlm(filename, comments=true) # e.g. "data/planck_Cl_TT_lowl.txt"
    l, Dl, ΔDlm, ΔDlp = data[:,1], data[:,2], data[:,3], data[:,4]
    return l, Dl, (ΔDlm, ΔDlp)
end

if false
    k = 10 .^ range(-4, 1, length=100) / Mpc # TODO: h in units
    plot(xlabel=L"\log_{10} \Big[ k / (h/\textrm{Mpc}) \Big]", ylabel=L"\log_{10} \Big[ P(k) / (\textrm{Mpc}/h)^3 \Big]")
    plot!(log10.(k / (par.h0/Mpc)), log10.(P.(Perturbations.(rec,k),0,k) / (Mpc/par.h0)^3)) # TODO
    vline!([log10(krm / (par.h0/Mpc))], color=:gray, linestyle=:dash)
    savefig("plots/power_spectrum_matter.pdf")
end

# Test source function
if false
    plot(xlabel=L"x = \log a", ylabel=L"S(x,k)", xlims=(-8,0))

    for k in [340*par.H0/c]
        perturb = Perturbations(rec, k)
        x = range(-20, 0, step=2*π*aH(par,0) / (20*c*k))
        plot!(x, S.(perturb, x) .* Cosmology.jl.(100, c*k*(η(bg,0.0) .- η.(bg,x))) / 1e-3, linewidth=0.5)
    end

    savefig("plots/source.pdf")
end

if false
    plot_Dl_against_Planck(:TT)
    plot_Dl_against_Planck(:TE)
    plot_Dl_against_Planck(:EE)
end

if true
    plot_Dl_varying_parameter(:h,           [0.57, 0.67, 0.77];       labelfunc = h    -> L"h = %$(h)")
    plot_Dl_varying_parameter(:Ωb0,         [0.02, 0.05, 0.08];       labelfunc = Ωb0  -> L"\Omega_{b0} = %$(Ωb0)")
    plot_Dl_varying_parameter(:Ωc0,         [0.217, 0.267, 0.317];    labelfunc = Ωc0  -> L"\Omega_{c0} = %$(Ωc0)")
    plot_Dl_varying_parameter(:Tγ0,         [2.6255, 2.7255, 2.8255]; labelfunc = Tγ0  -> L"T_{\gamma 0} = %$(Tγ0) \textrm{ K}")
    plot_Dl_varying_parameter(:Neff,        [0, 3.046, 6.092];        labelfunc = Neff -> L"N_\textrm{eff} = %$(Neff)")
    plot_Dl_varying_parameter(:ns,          [0.76, 0.96, 1.16];       labelfunc = ns   -> L"n_s = %$(ns)")
    plot_Dl_varying_parameter(:As,          [2e-8, 2e-9, 2e-10];      labelfunc = As   -> L"A_s = 2 \cdot 10^{%$(Int(round(log10(As/2))))}")
    plot_Dl_varying_parameter(:Yp,          [0, 0.24, 0.48];          labelfunc = Yp   -> L"Y_p = %$(Yp)")
    plot_Dl_varying_parameter(:z_reion_H,   [NaN, 8];                 labelfunc = z    -> isnan(z) ? "reionizatioff" : "reionization")
end

end
