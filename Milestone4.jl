module Milestone4

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM()
xrm = equality_rm(co)
krm = 1 / (c*η(co,xrm))

if false
    k = 10 .^ range(-4, 1, length=100) / Mpc # TODO: h in units
    plot(xlabel=L"\log_{10} \Big[ k / (h/\textrm{Mpc}) \Big]", ylabel=L"\log_{10} \Big[ P(k) / (\textrm{Mpc}/h)^3 \Big]")
    plot!(log10.(k / (co.h0/Mpc)), log10.(P.(co,0,k) / (Mpc/co.h0)^3)) # TODO
    vline!([log10(krm / (co.h0/Mpc))], color=:gray, linestyle=:dash)
    savefig("plots/power_spectrum_matter.pdf")
end

# Test source function
if true
    plot(xlabel=L"x = \log a", ylabel=L"S(x,k)", xlims=(-2,0))

    for k in [340*co.H0/c]
        x = range(-20, 0, step=2*π*aH(co,0) / (20*c*k))
        plot!(x, S.(co, x, k) .* Cosmology.jl.(100, c*k*(η(co,0.0) .- η.(co,x))) / 1e-3, linewidth=0.5)
    end

    savefig("plots/source.pdf")
end

# test Θl0
if true
    l = unique(Int.(round.(10 .^ range(0, 3, length=50))))
    plot(xlabel=L"\log_{10} l", ylabel=L"\log_{10} \Big[ \frac{l(l+1)}{2\pi} C_l \Big]")
    plot!(log10.(l), @. log10(l*(l+1)/(2*π)*Cl(co,l)))
    savefig("plots/power_spectrum_cmb.pdf")
end

end
