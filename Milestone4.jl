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

function lgrid(n1=10, n2=20, n3=150)
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

# test Θl0
if true
    plot_Dl_against_Planck(:TT)
    plot_Dl_against_Planck(:TE)
    plot_Dl_against_Planck(:EE)
end

end
