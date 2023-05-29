module Milestone3

using ForwardDiff

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants

par = Parameters()
bg = Background(par)
rec = Recombination(bg)
xrm = x_equality_rm(par)
xmΛ = x_equality_mΛ(par)
xrec = x_recombination(rec)

if true
    # analytical solutions from Dodelson eqs. (8.22), (8.31), (8.46), (8.52), (8.64)
    yf(x) = par.Ωm0/par.Ωr0 * a(x)
    Φy_small_anal(Φ0,y) = Φ0 / (10*y^3) * (16*√(1+y) + 9*y^3 + 2*y^2 - 8*y - 16)
    Φ_small_anal(p,x) = Φy_small_anal(Φ(p,-15.0), yf(x))
    xf(x,k) = c*k*η(bg,x) / √(3)
    Φ_large_anal(p,x) = p.k*Mpc > 0.05 ? 3 * (sin(xf(x,p.k))-xf(x,p.k)*cos(xf(x,p.k))) / xf(x,p.k)^3 * Φ(p,-15.0) : NaN

    δm_anal_matter(p,x) = x > xrm && c*p.k/aH(par,x) > 1 ? Φ(p,x) * 2*p.k^2*c^2*a(x) / (3*par.Ωm0*par.H0^2) : NaN
    A, B = 6.4, 0.44
    δc_anal_radiation(p,x) = x < xrm && B*c*p.k*η(bg,x) > 1 ? A * 3/2*Φ(p,-15.0) * (log(B*c*p.k*η(bg,x))) : NaN

    R(p,x) = 3/2 * Φ(p,x)
    Dm(y) = (y+2/3) * log((√(1+y)+1) / (√(1+y)-1)) - 2*√(1+y)
    Dp(y) = y + 2/3
    Dm′(y) = ForwardDiff.derivative(Dm, y)
    Dp′(y) = ForwardDiff.derivative(Dp, y)

    function δc_anal_Meszaros(p,x)
        a_H = a(x_horizon_entry(bg, p.k))
        y_H = a_H / a(xrm)
        y_m = 3*y_H
        if y_H > 1 || yf(x) < y_H || B*c*p.k*η(bg,x) < 1
            return NaN
        end
        C1, C2 = [Dp(y_m) Dm(y_m); Dp′(y_m) Dm′(y_m)] \ [A*R(p,-15.0)*log(B*y_m/y_H), A*R(p,-15.0)/y_m]
        return C1 * Dp(yf(x)) + C2 * Dm(yf(x))
    end

    fD(y) = 5*√(1+1/y) - 20/3*(1+1/y)^(3/2) + 8/3*((1+1/y)^(5/2)-1/y^(5/2))
    Θ0_plus_Φ_anal(p,x) = p.k*Mpc > 0.5 ? cos(p.k*c*η(bg,x)/√(3)) * exp(-p.k^2 * 3.1e6 * Mpc^2 * a(x)^(5/2) * fD(a(x)/a(xrm)) * (par.Ωb0*par.h0^2)^(-1) * (1-par.Yp/2)^(-1) * (par.Ωm0*par.h0^2)^(-1/2)) * (Φ(p,-10.7)+Θl(p,-10.7,0)) : NaN

    series = [
        ("plots/ThetaP.pdf", Dict(:ylabel => L"\Theta^P_l"), [
            ((p,x) -> ΘPl(p,x,0), Dict(:linestyle => :solid, :label => L"l=0")),
            ((p,x) -> ΘPl(p,x,1), Dict(:linestyle => :dash, :label => L"l=1")),
            ((p,x) -> ΘPl(p,x,2), Dict(:linestyle => :dot, :label => L"l=2"))
        ]),

        ("plots/oscillations.pdf", Dict(:ylabel => L"\Theta_0+\Phi"), [
            ((p,x) -> Θl(p,x,0) + Φ(p,x), Dict(:linestyle => :solid, :label => nothing)),
            ((p,x) -> Θ0_plus_Φ_anal(p,x), Dict(:linestyle => :solid, :color => :gray, :alpha => 0.5, :label => nothing)),
        ]),

        ("plots/overdensity.pdf", Dict(:ylabel => L"\log_{10}|\delta|"), [
            ((p,x) -> log10(abs(δm_anal_matter(p,x))), Dict(:linestyle => :solid, :color => :gray, :alpha => 0.5, :label => nothing)),
            ((p,x) -> log10(abs(δc_anal_radiation(p,x))), Dict(:linestyle => :solid, :color => :gray, :alpha => 0.5, :label => nothing)),
            ((p,x) -> log10(abs(δc_anal_Meszaros(p,x))), Dict(:linestyle => :solid, :color => :gray, :alpha => 0.5, :label => nothing)),
            ((p,x) -> log10(abs(δc(p,x))), Dict(:linestyle => :solid, :label => L"\delta=\delta_c")),
            ((p,x) -> log10(abs(δb(p,x))), Dict(:linestyle => :dash,  :label => L"\delta=\delta_b")),
            ((p,x) -> log10(abs(4*Θl(p,x,0))), Dict(:linestyle => :dashdot,  :label => L"\delta=\delta_\gamma=4\Theta_0")),
            ((p,x) -> log10(abs(4*Nl(p,x,0))), Dict(:linestyle => :dot,      :label => L"\delta=\delta_\nu=4\mathcal{N}_0")),
        ]),

        ("plots/potentials.pdf", Dict(:ylabel => L"\{\Phi,\Psi\}", :ylims => (-0.8, +0.8)), [
            ((p,x) -> Φ_small_anal(p,x), Dict(:linestyle => :solid, :color => :gray, :alpha => 0.5, :label => nothing)),
            ((p,x) -> Φ_large_anal(p,x), Dict(:linestyle => :solid, :color => :gray, :alpha => 0.5, :label => nothing)),
            ((p,x) -> Φ(p,x), Dict(:linestyle => :dash, :label => L"\Phi")),
            ((p,x) -> Ψ(p,x), Dict(:linestyle => :dot, :label => L"\Psi")),
            ((p,x) -> Φ(p,x) + Ψ(p,x), Dict(:linestyle => :solid, :label => L"\Phi+\Psi"))
        ]),

        ("plots/velocity.pdf", Dict(:ylabel => L"\log_{10}|v|"), [
            ((p,x) -> log10(abs(vc(p,x))), Dict(:linestyle => :solid, :label => L"v=v_c")),
            ((p,x) -> log10(abs(vb(p,x))), Dict(:linestyle => :dash, :label => L"v=v_b")),
            ((p,x) -> log10(abs(-3*Θl(p,x,1))), Dict(:linestyle => :dashdot, :label => L"v=v_\gamma=-3\Theta_1")),
            ((p,x) -> log10(abs(-3*Nl(p,x,1))), Dict(:linestyle => :dot,     :label => L"v=v_\gamma=-3\mathcal{N}_1")),
        ]),

        ("plots/ThetalN0.pdf", Dict(:ylabel => L"\{\Theta_0,\mathcal{N}_0\}"), [((p,x) -> Θl(p,x,0), Dict(:linestyle => :solid, :label => L"\Theta_0")), ((p,x) -> Nl(p,x,0), Dict(:linestyle => :dash, :label => L"\mathcal{N}_0"))]),
        ("plots/ThetalN1.pdf", Dict(:ylabel => L"\{\Theta_1,\mathcal{N}_1\}"), [((p,x) -> Θl(p,x,1), Dict(:linestyle => :solid, :label => L"\Theta_1")), ((p,x) -> Nl(p,x,1), Dict(:linestyle => :dash, :label => L"\mathcal{N}_1"))]),
        ("plots/ThetalN2.pdf", Dict(:ylabel => L"\{\Theta_2,\mathcal{N}_2\}"), [((p,x) -> Θl(p,x,2), Dict(:linestyle => :solid, :label => L"\Theta_2")), ((p,x) -> Nl(p,x,2), Dict(:linestyle => :dash, :label => L"\mathcal{N}_2"))]),
    ]

    for (ploti, (filename, plotsettings, func_linesettings)) in enumerate(series)
        println("Plotting $filename")
        plot(xlabel=L"x = \log a", xlims=(-15, 0), xticks=-25:1:5, legend_position=:bottomleft; plotsettings...)
        for (i, k) in enumerate([1e-0, 1e-1, 1e-2, 1e-3] / Mpc)
            perturb = PerturbationMode(rec, k)
            xhor = x_horizon_entry(bg, k)
            x = Cosmology.extend(Cosmology.integration_points(perturb.y), 3)
            for (func, linesettings) in func_linesettings
                plot!(x, func.(perturb, x), alpha=1.0, linewidth=0.5, color=i; linesettings..., label=nothing)
            end
            if ploti == 1
                # put k-labels and k-independent time marks on first plot only
                vline!([-21], color=i, label=L"k=10^{%$(Int(round(log10(k*Mpc))))}/\textrm{Mpc}") # label each k-value once
            end
            vline!([xhor], color=i, linestyle=:solid; alpha=0.3, linewidth=3.0, label=nothing)
            vline!([xrm, xmΛ, xrec], color=:gray, linestyle=:solid; alpha=0.3, linewidth=3.0, label=nothing)
        end

        # add quantity (if more than one so ambiguous) to legend with a black dummy plot
        for (func, linesettings) in func_linesettings
            vline!([-21], color=:black; linesettings...) # dummy outside plot area
        end
        savefig(filename)
    end
end

end
