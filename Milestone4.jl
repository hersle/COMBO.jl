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

Dl(l, Cl) = l * (l+1) / (2*π) * Cl
logerrm(y, Δy) = -log10(1-Δy/y)
logerrp(y, Δy) = +log10(1+Δy/y)

function plot_Cl(ls, Cls, datafilename, polarization, filename)
    plot(xlabel=L"\log_{10} l", ylabel=L"\log_{10} \Big[ \frac{l(l+1)}{2\pi} C_l^{%$(polarization)} \Big]")

    plot!(log10.(ls), log10.(Dl.(ls, Cls)), label="ΛCDM prediction")

    ls_data, Dls_data, ΔDlms_data, ΔDlps_data = read_planck_data(datafilename)
    scatter!(log10.(ls_data), log10.(Dls_data); yerror=(logerrm.(Dls_data, ΔDlms_data), logerrp.(Dls_data, ΔDlps_data)), color=:black, markersize=1, label="Planck measurements")
    # TODO: don't do log units?

    savefig(filename)
end

function read_planck_data(filename)
    data = readdlm(filename, comments=true) # e.g. "data/planck_Cl_TT_lowl.txt"
    l, Dl, ΔDlm, ΔDlp = data[:,1], data[:,2], data[:,3], data[:,4]
    Dl   .*= (1e-6 / par.Tγ0)^2 # convert to my preferred units
    ΔDlm .*= (1e-6 / par.Tγ0)^2
    ΔDlp .*= (1e-6 / par.Tγ0)^2
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
if false
    l = unique(Int.(round.(10 .^ range(0, 3.4, length=300))))
    Cl = Cls(rec,l)

    plot_Cl(ls, Cls, "plots/power_spectrum_cmb.pdf")
end

end
