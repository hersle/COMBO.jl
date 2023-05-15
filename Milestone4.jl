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

Dl(l, Cl) = l * (l+1) / (2*π) * Cl * (par.Tγ0 / 1e-6)^2 # convert to "Planck units"
logerrm(y, Δy) = -log10(1-Δy/y)
logerrp(y, Δy) = +log10(1+Δy/y)

function plot_Cl(ls, Cls, datafilename, polarization, filename)
    plot(xlabel=L"\log_{10} l", ylabel=L"C_l^\mathrm{%$(polarization)} \, T_{\gamma 0} \cdot \frac{l(l+1)}{2\pi} \,/\, (\mathrm{\mu} K)^2", legend_position=:topleft)

    plot!(log10.(ls), Dl.(ls, Cls), label="ΛCDM prediction", marker=:circle, markersize=1.0, markerstrokecolor=1)

    if !isnothing(datafilename)
        ls_data, Dls_data, ΔDlms_data, ΔDlps_data = read_planck_data(datafilename)
        scatter!(log10.(ls_data), Dls_data; yerror=(ΔDlms_data, ΔDlps_data), color=:black, marker=:square, markersize=1, label="Planck measurements")
    end

    savefig(filename)
end

function read_planck_data(filename)
    data = readdlm(filename, comments=true) # e.g. "data/planck_Cl_TT_lowl.txt"
    l, Dl, ΔDlm, ΔDlp = data[:,1], data[:,2], data[:,3], data[:,4]
    #Dl   .*= (1e-6 / par.Tγ0)^2 # for converting to "dimensionless units"
    #ΔDlm .*= (1e-6 / par.Tγ0)^2
    #ΔDlp .*= (1e-6 / par.Tγ0)^2
    return l, Dl, ΔDlm, ΔDlp
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

# test Θl0
if true
    l = unique(Int.(round.(10 .^ range(0, 3.4, length=200))))
    plot_Cl(l, Cl(rec,l,:TT), "data/planck_Cl_TT.txt", "TT", "plots/power_spectrum_CMB_TT.pdf")
    plot_Cl(l, Cl(rec,l,:TE), "data/planck_Cl_TE.txt", "TE", "plots/power_spectrum_CMB_TE.pdf")
    plot_Cl(l, Cl(rec,l,:EE), "data/planck_Cl_EE.txt", "EE", "plots/power_spectrum_CMB_EE.pdf")
end

end
