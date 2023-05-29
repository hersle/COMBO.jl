module Milestone4

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using DelimitedFiles

par = Parameters()
bg = Background(par)
rec = Recombination(bg)
xrm = x_equality_rm(par)

function plot_Dl(series, polarization, filename)
    plot(xlabel=L"\log_{10} l", ylabel=L"D_l^\mathrm{%$(polarization)} / (\mathrm{\mu K})^2", xlims=(0, 3.4), legend_position=:topleft)

    for (i, (l, Dl, ΔDl, settings)) in enumerate(series)
        plot!(log10.(l), Dl; yerror=ΔDl, color=i, markerstrokecolor=i, markersize=1, marker=:circle, settings...)
    end

    savefig(filename)
end

function plot_Dl_against_Planck(type)
    pspec = CMBPowerSpectrum(rec, type; n1=10, n2=20, n3=150)
    plot_Dl([
        (pspec.l, pspec.Dl/1e-6^2, nothing, Dict(:label => "Our ΛCDM prediction")),
        (read_Planck_data("data/Planck_Dl_$type.txt")..., Dict(:label => "Planck's measurements", :color => :black, :marker => :square, :markerstrokecolor => :black, :markersize => 0.75, :markerstrokewidth => 0.50, :seriestype => :scatter)),
    ], "$type", "plots/power_spectrum_CMB_$type.pdf")
end

function plot_Dl_varying_parameter(paramkey, values; labelfunc=val -> nothing, type=:TT)
    series = []
    for (i, value) in enumerate(values)
        rec = Recombination(Background(Parameters(; Dict(paramkey => value)...)))

        # * With 3 values, the idea is always to plot
        #   the Planck values[2] in "neutral" black,
        #   the lower values[1] in "colder" blue,
        #   the higher values[3] in "warmer" red
        # * With 2 values, the idea is to plot
        #   a   feature  enabled with values[2] in "neutral" black,
        #   the feature disabled with values[1] in "cold" blue
        color = [:blue, :black, :red][i]

        pspec = CMBPowerSpectrum(rec, type; n1=10, n2=20, n3=150)
        push!(series, (pspec.l, pspec.Dl/1e-6^2, nothing, Dict(:color => color, :marker => :none, :label => labelfunc(value))))
    end
    plot_Dl(series, "$type", "plots/power_spectrum_CMB_varying_$paramkey.pdf")
end

function read_Planck_data(filename)
    data = readdlm(filename, Float64, comments=true) # e.g. "data/planck_Cl_TT_lowl.txt"
    l, Dl, ΔDlm, ΔDlp = data[:,1], data[:,2], data[:,3], data[:,4]
    return l, Dl, (ΔDlm, ΔDlp)
end

function read_Pk_data(filename; last_column_is_upper_bound=false)
    data = readdlm(filename, comments=true) # e.g. "data/planck_Cl_TT_lowl.txt"
    k, Pk, lastcol = data[:,1], data[:,2], data[:,3]
    ΔPk = last_column_is_upper_bound ? lastcol .- Pk : lastcol
    return k, Pk, ΔPk
end

if true
    k = 10 .^ range(-4, +2.0, length=200) / Mpc
    pspec = MatterPowerSpectrum(rec, k)

    plot(xlabel=L"\log_{10} \Big[ k / (h/\textrm{Mpc}) \Big]", ylabel=L"\log_{10} \Big[ P(k) / (\textrm{Mpc}/h)^3 \Big]", xlims=(-3, 1), ylims=(0, 5), legend_position=:topright)
    plot!(log10.(pspec.k/(par.h0/Mpc)), log10.(pspec.Pk/(Mpc/par.h0)^3), label="Our ΛCDM prediction")

    dataseries = [
        ("data/Pk_SDSS_DR7_LRG.txt", false, "SDSS DR7 LRG data"),
        ("data/Pk_WMAP_ACT.txt", true, "WMAP & ACT data"),
        ("data/Pk_Lyman_alpha.txt", true, "eBOSS data")
    ]
    for (i, (filename, last_column_is_upper_bound, label)) in enumerate(dataseries)
        k, Pk, ΔPk = read_Pk_data(filename; last_column_is_upper_bound=last_column_is_upper_bound)
        scatter!(log10.(k), log10.(Pk), yerror=(log10.(Pk) .- log10.(Pk-ΔPk), log10.(Pk+ΔPk) .- log10.(Pk)), marker=:square, markersize=1.50, markerstrokewidth=1.0, color=i+1, markerstrokecolor=i+1, label=label)
    end

    vline!(log10.([1 / (c*η(bg,xrm)), a(xrm)*H(par,xrm)/c] ./ (par.h0/Mpc)), color=:gray, linestyle=[:dash, :dot], label=[L"k=1/c \, \eta_\textrm{eq}", L"k=\mathcal{H}(a_\textrm{eq}) / c"])
    savefig("plots/power_spectrum_matter.pdf")
end

if true
    @time plot_Dl_against_Planck(:TT)
    @time plot_Dl_against_Planck(:TE)
    @time plot_Dl_against_Planck(:EE)
end

if true
    plot_Dl_varying_parameter(:h,           [0.57, 0.67, 0.77];        labelfunc = h    -> L"h = %$(h)")
    plot_Dl_varying_parameter(:Ωb0,         [0.02, 0.05, 0.08];        labelfunc = Ωb0  -> L"\Omega_{b0} = %$(Ωb0)")
    plot_Dl_varying_parameter(:Ωc0,         [0.217, 0.267, 0.317];     labelfunc = Ωc0  -> L"\Omega_{c0} = %$(Ωc0)")
    plot_Dl_varying_parameter(:Tγ0,         [2.6255, 2.7255, 2.8255];  labelfunc = Tγ0  -> L"T_{\gamma 0} = %$(Tγ0) \textrm{ K}")
    plot_Dl_varying_parameter(:Neff,        [0, 3.046, 6.092];         labelfunc = Neff -> L"N_\textrm{eff} = %$(Neff)")
    plot_Dl_varying_parameter(:ns,          [0.765, 0.965, 1.165];     labelfunc = ns   -> L"n_s = %$(ns)")
    plot_Dl_varying_parameter(:As,          [1.05e-9, 2.1e-9, 4.2e-9]; labelfunc = As   -> L"A_s = %$(As/1e-9) \cdot 10^{-9}")
    plot_Dl_varying_parameter(:Yp,          [0, 0.245, 0.49];          labelfunc = Yp   -> L"Y_p = %$(Yp)")
    plot_Dl_varying_parameter(:z_reion_H,   [NaN, 8];                  labelfunc = z    -> isnan(z) ? "reionizatioff" : "reionization")
    #plot_Dl_varying_parameter(:Ωk0,         [-0.5, 0, +0.5];           labelfunc = Ωk0  -> L"\Omega_{k0} = %$(Ωk0)") # TODO: need to implement curvature in perturbations first
end

if true
    ls = [1, 10, 100, 1000]
    η0 = η(bg,0)

    # plot dΘl0_dl = S * j?
    plot(xlabel=L"x = \log a", ylabel=L"\partial \Theta_{10}(x,k) / \partial x", legend_position=:topright)
    kcη0s = [10, 100, 1000]
    xs = range(-10, 0, step=0.01) # TODO: 0.02?
    ks = kcη0s / (c*η0)
    STs = Cosmology.grid_S(rec, Cosmology.ST, xs, ks; spline_first=false)
    χs = χ.(bg, xs) # = c * (η(0) - η(x))
    dΘTl0_dxs = dΘTl0_dx(10, xs, ks, STs, χs)
    ys = [dΘTl0_dxs[:,i] for i in 1:length(ks)]
    plot!(xs, ys, linewidth=0.3, xlims=(-8, 0), ylims=(-1, +1), yticks=-1.0:0.20:+1.0, label=reshape([L"k=%$(kcη0)/(c\eta_0)^{-1}" for kcη0 in kcη0s], 1, :))
    plot!(xs, ys, linewidth=0.3, xlims=(-1, 0), ylims=(-0.03, +0.03), yticks=-0.03:0.01:+0.03, xticks=-1:0.1:0, subplot=2, inset=(1, bbox(0.10, 0.43, 0.7, 0.5, :right)), label=nothing)
    savefig("plots/dThetal0_dx.pdf")


    # plot Θl0 # TODO: ΘEl0, too?
    plot(xlabel=L"k / (c \eta_0)^{-1}", ylabel=L"\Theta_l(x=0,k)", xticks=0:100:4000, ylims=(-0.1, +0.1), legend_position=:topright)
    ks = range(1/(c*η0), 1200/(c*η0), step=2*π/(c*η0*10))
    χs = χ.(bg, xs) # = c * (η(0) - η(x))
    STs = Cosmology.grid_S(rec, Cosmology.ST, xs, ks; spline_first=true)
    for l in ls
        plot!(ks*c*η0, ΘTl0(l, xs, ks, STs, χs), linewidth=0.3, label=L"l=%$l")
    end
    savefig("plots/Theta0.pdf")


    # plot dCl_dk
    for (i, l) in enumerate(ls)
        kcη0min = max(0, l-100)
        kcη0max = kcη0min + 300
        y = dCl_dk_TT(l,xs,ks,STs,nothing,χs,par) / (c*η0)
        e10 = Int(round(log10(maximum(y))))
        y ./= 10.0^e10
        plot(xlabel=L"k / (c \eta_0)^{-1}", ylabel=L"\mathrm{d} C_l(k) / \mathrm{d} k / 10^{%$(e10)} (c \eta_0)^{-1}", xlims=(kcη0min, kcη0max), xticks=0:100:2000, legend_position=:topright)
        plot!(c*ks*η0, y, linewidth=0.5, color=i, label=L"l=%$l")
        savefig("plots/dCl_dk_$l.pdf")
    end
end

end
