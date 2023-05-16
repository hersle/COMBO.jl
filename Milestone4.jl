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
    @time plot_Dl_against_Planck(:TT)
    @time plot_Dl_against_Planck(:TE)
    @time plot_Dl_against_Planck(:EE)
end

if false
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

if true
    ls = 1, 10, 100, 1000
    η0 = η(bg,0)
    x = unique(vcat(range(-8, -3, length=1000), range(-3, 0, length=4000)))
    k = range(1/(c*η0), 1200/(c*η0), length=3000)
    Sspl, SEspl = Cosmology.spline_S(rec, [Cosmology.S, Cosmology.SE])
    Θl0(l,k) = Cosmology.Θl0(l, k, Sspl, bg.η)
    dCl_dk(l,k) = Cosmology.dCl_dk_TT(l,k,Sspl,SEspl,bg.η,par)

    # plot dΘl0_dl = S * j?
    plot(xlabel=L"x = \log a", ylabel=L"\partial \Theta_{l}(x,k) / \partial x", legend_position=:bottomright)
    kcη0 = [400, 4000]
    ks = kcη0 / (c*η0)
    kMpc = ks / Mpc
    ys = [Cosmology.dΘl0_dx.(10, k, x, Sspl, Ref(bg.η)) for k in ks]
    plot!(x, ys, linewidth=0.3, xlims=(-8, 0), ylims=(-0.02, +0.02), yticks=-0.02:0.01:+0.02, label=[L"k=400/(c \eta_0)"  L"k=4000/(c \eta_0)"])
    plot!(x, ys, linewidth=0.3, xlims=(-1, 0), ylims=(-0.02, +0.02), yticks=-0.02:0.01:+0.02, xticks=-1:0.1:0, subplot=2, inset=(1, bbox(0.15, 0.03, 0.7, 0.5, :right)), label=nothing)
    savefig("plots/dThetal0_dx.pdf")

    # plot Θl0 # TODO: ΘEl0?
    # TODO: integrate Θl0 from k=0, forcing it to start from 0?
    plot(xlabel=L"k / c \eta_0", ylabel=L"\Theta_l(x=0,k)", xticks=0:100:4000, ylims=(-0.1, +0.1), legend_position=:topright)
    for l in ls
        plot!(k*c*η0, Θl0.(l, k), linewidth=0.3, label=L"l=%$l")
    end
    savefig("plots/Theta0.pdf")

    k = range(1/(c*η0), 300/(c*η0), length=2000)
    for (i, l) in enumerate(ls)
        plot(xlabel=L"k / c \eta_0", ylabel=L"\mathrm{d} C_l(k) / \mathrm{d} k / (c \eta_0)^{-1}", xlims=(0, 300), xticks=0:100:500, legend_position=:topright)
        plot!(k*c*η0, dCl_dk.(l,k) / (c*η0), linewidth=0.5, color=i, label=L"l=%$l")
        savefig("plots/dCl_dk_$l.pdf")
    end

    # plot dCl_dk
end

end
