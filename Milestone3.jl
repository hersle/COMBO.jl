module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings
using DifferentialEquations

#co = ΛCDM(h=0.7, Neff=0, Ωb0=0.05, Ωc0=0.45, Yp=0, z_reion_H=NaN) # to compare with Hans
#println("Tight coupling end (k = $(k*Mpc)/Mpc): ", format_time_variations(co, xtce))

#=
if false
    Cosmology.perturbations_splines(co)
    Cosmology.perturbations_mode(co, 1 / Mpc; tight=true, solver=KenCarp4(autodiff=false))
    Cosmology.perturbations_mode(co, 1 / Mpc; tight=true, solver=KenCarp4(autodiff=false))
    Cosmology.perturbations_mode(co, 1 / Mpc; tight=true, solver=Tsit5())
    Cosmology.perturbations_mode(co, 1 / Mpc; tight=true, solver=Tsit5())
end
=#

# Test splining
#=
if false
    kmin, kmax = 0.00005 / Mpc, 1.0 / Mpc
    xmin, xmax = -20.0, 0.0

    series = [
        ("logarithmic", 10 .^ (log10(kmin) .+ (log10(kmax)-log10(kmin)) * range(0, 1; length=200))),
        ("quartic", kmin .+ (kmax-kmin) * range(0, 1; length=200) .^ 4),
        ("cubic", kmin .+ (kmax-kmin) * range(0, 1; length=200) .^ 3),
        ("quadratic", kmin .+ (kmax-kmin) * range(0, 1; length=200) .^ 2),
        ("linear", kmin .+ (kmax-kmin) * range(0, 1; length=200) .^ 1),
    ]

    for (variant, splineks) in series
        filename = "plots/splinetest_$variant.pdf"
        println("Plotting $filename")

        ks = 10 .^ (log10(kmin) .+ (log10(kmax)-log10(kmin)) * sort(rand(100)))
        xs = xmin .+ (xmax-xmin) * sort(rand(341)) .^ 1 # linear spacing

        y1s = Cosmology.perturbations_splines(co; ks=splineks)

        maxdy = zeros(length(xs), length(ks))
        for (j, k) in enumerate(ks)
            y1spls = [x -> spl(x, k) for spl in y1s] # list of functions
            _, y2spls = Cosmology.perturbations_mode(co, k) # list of functions
            for l in 1:length(y1spls)
                y1 = y1spls[l].(xs)
                y2 = y2spls[l].(xs)
                y1 /= maximum(abs.(y1))
                y2 /= maximum(abs.(y2)) # TODO: divide by max(y1) or max(y2)?
                for (i, x) in enumerate(xs)
                    maxdy[i,j] = max(maxdy[i,j], abs(y1[i]-y2[i]))
                end
            end
        end

        x = xs
        y = log10.(ks * Mpc)
        z = log10.(maxdy') # transpose
        plot(title=L"\log_{10} \left( \max_i |\hat{y}_i(x,k)-\hat{y}^k_i(x)| \right) \quad \textrm{(%$(variant) spacing)}", xlabel=L"x=\log a", ylabel=L"\log_{10} \left( k \cdot \textrm{Mpc} \right)", xlims=(x[1], x[end]), ylims=(y[1], y[end]), grid=nothing, colorbar_ticks=-20:0, clims=(-8, 0))
        heatmap!(x, y, z, c=cgrad(:Reds), label=L"\log_{10} dy")
        savefig(filename)
    end
end
=#

if true
    par = Parameters()
    bg = Background(par)
    rec = Recombination(bg)
    xrm = equality_rm(par)
    xmΛ = equality_mΛ(par)
    xrec = time_recombination(rec)

    # analytical solutions
    yf(x) = par.Ωm0/par.Ωr0 * a(x)
    Φy_small_anal(Φ0,y) = Φ0 / (10*y^3) * (16*√(1+y) + 9*y^3 + 2*y^2 - 8*y - 16)
    dΦy_small_anal(Φ0,y; ϵ=1e-4) = (Φy_small_anal(Φ0,y+ϵ/2)-Φy_small_anal(Φ0,y-ϵ/2)) / ϵ
    Φ_small_anal(p,x) = Φy_small_anal(Φ(p,-15.0), yf(x))
    dΦ_small_anal(p,x) = dΦy_small_anal(Φ(p,-15.0), yf(x))
    xf(x,k) = c*k*η(bg,x) / √(3)
    Φ_large_anal(p,x) = p.k*Mpc > 0.05 ? 3 * (sin(xf(x,p.k))-xf(x,p.k)*cos(xf(x,p.k))) / xf(x,p.k)^3 * Φ(p,-15.0) : NaN

    δm_anal_matter(p,x) = x > xrm && c*p.k/aH(par,x) > 1 ? Φ(p,x) * 2*p.k^2*c^2*a(x) / (3*par.Ωm0*par.H0^2) : NaN
    A, B = 6.4, 0.44
    δc_anal_radiation(p,x) = x < xrm && B*c*p.k*η(bg,x) > 1 ? A * 3/2*Φ(p,-15.0) * (log(B*c*p.k*η(bg,x))) : NaN
    δc_anal_small(p,x) = x < xrm ? 6*(yf(x)+1)/(3*yf(x)+4) * (yf(x)*dΦ_small_anal(p,x)+Φ_small_anal(p,x)) : NaN # Dodelson (8.22)

    R(p,x) = 3/2 * Φ(p,x)
    Dm(y) = (y+2/3) * log((√(1+y)+1) / (√(1+y)-1)) - 2*√(1+y)
    Dp(y) = y + 2/3
    dDm(y; ϵ=1e-4) = (Dm(y+ϵ/2) - Dm(y-ϵ/2)) / ϵ # derivatives; forward-difference needed for matching
    dDp(y; ϵ=1e-4) = (Dp(y+ϵ/2) - Dp(y-ϵ/2)) / ϵ # derivatives; forward-difference needed for matching

    function δc_anal_Meszaros(p,x)
        a_H = a(time_horizon_entry(bg, p.k))
        y_H = a_H / a(xrm)
        y_m = 3*y_H
        if y_H > 1 || yf(x) < y_H || B*c*p.k*η(bg,x) < 1
            return NaN
        end
        C1, C2 = [Dp(y_m) Dm(y_m); dDp(y_m) dDm(y_m)] \ [A*R(p,-15.0)*log(B*y_m/y_H), A*R(p,-15.0)/y_m]
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
            perturb = Perturbations(rec, k)
            xhor = time_horizon_entry(bg, k)
            x = Cosmology.extendx(Cosmology.splinex(perturb.qty_splines[1]), 3)
            for (func, linesettings) in func_linesettings
                plot!(x, func.(perturb, x), alpha=1.0, linewidth=0.5, color=i; linesettings..., label=nothing)
            end
            if ploti == 1
                # put k-labels and k-independent time marks on first plot only
                vline!([-21], color=i, label=L"k=10^{%$(Int(round(log10(k*Mpc))))}/\textrm{Mpc}") # label each k-value once
                vline!([xrm, xmΛ], color=:gray, linestyle=:solid; alpha=0.3, linewidth=3.0, label=nothing)
                vline!([xrec], color=:gray, linestyle=:solid; alpha=0.3, linewidth=3.0, label=nothing)
            end
            #vline!([time_tight_coupling(co, k)], color=:gray, linestyle=:dash; linewidth=0.5, label=nothing)
            vline!([xhor], color=i, linestyle=:solid; alpha=0.3, linewidth=3.0, label=nothing)
        end

        # add quantity (if more than one so ambiguous) to legend with a black dummy plot
        for (func, linesettings) in func_linesettings
            vline!([-21], color=:black; linesettings...) # dummy outside plot area
        end
        savefig(filename)
    end
end

#=
if true
    # plot Θ0 zoomed in (to check whether we capture the most rapid oscillations)
    filename = "plots/Thetal0_zoom.pdf"
    println("Plotting $filename")
    k = ks[2]
    _, y1 = Cosmology.perturbations_mode(rec, k; tight=true ) # 1D (x) splines for each "exact" k for tight+full system
    x, y2 = Cosmology.perturbations_mode(rec, k; tight=false) # 1D (x) splines for each "exact" k for       full system
    x = x[x .> -1.1]
    x = Cosmology.extendx(x, 20)
    func(x, y) = y[Cosmology.i_Θl(0)](x)
    plot(xlabel=L"x = \log a", ylabel=L"\Theta_0", xlims=(-1.0, 0), xticks=-5.0:1:0, legend_position=:topright)
    plot!(x, func.(x, Ref(y1)), label=L"\texttt{ Tsit5}    \textrm{ method with tight coupling}")
    plot!(x, func.(x, Ref(y2)), label=L"\texttt{ KenCarp4} \textrm{ method without tight coupling}")
    savefig(filename)

end

# TODO: compare RK methods. plot max(yM1 - yM2) for each x and some k
if true
    ks = [1e-1, 1e-2, 1e-3] / Mpc
    filename = "plots/perturbation_methods.pdf"
    println("Plotting $filename")
    plot(xlabel=L"x = \log a", ylabel=L"\log_{10} \Big( \max_i{|\hat{y}^{1}_i-\hat{y}^{2}_i|} \Big)", xlims=(-15, 0), xticks=-25:5:5, legend_position=:topleft)
    for (i_k, k) in enumerate(ks)
        x, y1s  = Cosmology.perturbations_mode(co, k; tight=false, solver=KenCarp4(autodiff=false))
        _, y2s  = Cosmology.perturbations_mode(co, k; tight=true,  solver=KenCarp4(autodiff=false))
        _, y3s  = Cosmology.perturbations_mode(co, k; tight=true,  solver=Tsit5())
        x = Cosmology.extendx(x, 3) # plot with extra points between every spline point for more smoothness
        y1s = [y1.(x) / maximum(abs.(y1.(x))) for y1 in y1s] # [i_q, i_x]
        y2s = [y2.(x) / maximum(abs.(y2.(x))) for y2 in y2s] # [i_q, i_x]
        y3s = [y3.(x) / maximum(abs.(y3.(x))) for y3 in y3s] # [i_q, i_x]
        dy12 = [abs.(y1 .- y2) for (y1, y2) in zip(y1s, y2s)]
        dy13 = [abs.(y1 .- y3) for (y1, y3) in zip(y1s, y3s)]
        dy12 = [maximum([dy12[i][j] for i in 1:length(dy12)]) for j in 1:length(x)]
        dy13 = [maximum([dy13[i][j] for i in 1:length(dy13)]) for j in 1:length(x)]
        plot!(x, log10.(dy13); alpha=0.5, color=1+i_k, label=nothing)
        plot!(x, log10.(dy12); alpha=1.0, color=1+i_k, label=L"k=10^{%$(Int(round(log10(k*Mpc))))}/\textrm{Mpc}")
    end
    vline!([-21], color=:black, alpha=1.0, label=L"\textrm{\texttt{KenCarp4} (full) vs. \texttt{KenCarp4} (tight+full)}")
    vline!([-21], color=:black, alpha=0.5, label=L"\textrm{\texttt{KenCarp4} (full) vs. \texttt{Tsit5} (tight+full)}")
    savefig(filename)
end
=#

end
