using Profile
using BenchmarkTools

const polarization = true
const neutrinos = true

const lmax = 10 # use >= 6 without neutrinos/polarization; >= 10 with neutrinos/polarization
@assert lmax >= 4 # equations are ambiguous with lmax <= 3 (how to handle l<2, l=2 and l=lmax?)

# variable index map (1-based),
# ordered so full system *extends* the tightly coupled one (useful when switching)
const i_δc = 1
const i_vc = 2
const i_δb = 3
const i_vb = 4
const i_Φ  = 5
const i_Nl(l) = 6 + l
const i_Θl(l) = i_Nl(lmax) + 1 + l
const i_ΘPl(l) = i_Θl(lmax) + 1 + l # variables here and below are only part of the full system

const i_max_tight = i_Θl(1)     # last variable of tight system
const i_max_full  = i_ΘPl(lmax) # last variable of full system

struct PerturbationMode
    rec::Recombination
    k::Float64
    y::ODESolution

    function PerturbationMode(rec::Recombination, k::Float64; tight=false, kwargs...)
        new(rec, k, perturbations_mode(rec, k; tight=tight, kwargs...))
    end
end

Base.broadcastable(pert::PerturbationMode) = Ref(pert)

function time_tight_coupling(rec::Recombination, k; tol::Float64=10.0)
    x1 = find_zero(x -> abs(dτ(rec,x)) - tol,                        (-20.0, +20.0)) # tight coupling assumes     1/τ′ << 1
    x2 = find_zero(x -> abs(dτ(rec,x)) - tol * c*k/aH(rec.bg.par,x), (-20.0, +20.0)) # tight coupling assumes ck/̇(aτ′) << 1
    x3 = time_switch_Peebles(rec.bg.par) # switch no later than when recombination begins # TODO: or -8.3?
    x4 = -10.0 # switch before this to avoid kink in ΘP # TODO: find dynamic way of computing this?
    return min(x1, x2, x3, x4) # take the earliest of these times
end

function time_horizon_entry(bg::Background, k)
    return find_zero(x -> k * c*η(bg,x) - 1, (-20.0, +20.0))
end

function perturbations_initial_conditions(par::Parameters, η, τ, x0, k)
    fν = par.Ων0 / par.Ωr0
    Ψ  = -1 / (3/2 + 2*fν/5)
    τ′x = ForwardDiff.derivative(τ, x0)

    y = Vector{Float64}(undef, i_max_full)
    y[i_δc] = y[i_δb] = -3/2 * Ψ
    y[i_vc] = y[i_vb] = -k*c/(2*aH(par,x0)) * Ψ
    y[i_Φ]            = -(1 + 2*fν/5) * Ψ
    y[i_Θl(0)]        = -1/2 * Ψ
    y[i_Θl(1)]        = -c*k / (3*aH(par,x0)) * y[i_Θl(0)]
    y[i_Θl(2)]        = (polarization ? -8/15 : -20/45) * c*k / (aH(par,x0)*τ′x) * y[i_Θl(1)]
    for l in 3:lmax
        y[i_Θl(l)]    = -l/(2*l+1) * c*k/(aH(par,x0)*τ′x) * y[i_Θl(l-1)] # recursive relation
    end

    if polarization
        y[i_ΘPl(0)]       = 5/4 * y[i_Θl(2)]
        y[i_ΘPl(1)]       = -c*k/(4*aH(par,x0)*τ′x) * y[i_Θl(2)]
        y[i_ΘPl(2)]       = 1/4 * y[i_Θl(2)]
        for l in 3:lmax
            y[i_ΘPl(l)]    = -l/(2*l+1) * c*k/(aH(par,x0)*τ′x) * y[i_ΘPl(l-1)] # recursive relation
        end
    else
        y[i_ΘPl(0):i_ΘPl(lmax)] .= 0.0
    end

    if neutrinos
        y[i_Nl(0)] = -1/2 * Ψ
        y[i_Nl(1)] = +c*k/(6*aH(par,x0)) * Ψ
        y[i_Nl(2)] = +(c*k*a(x0)/par.H0)^2 / (30*par.Ωr0) * Ψ # expand to avoid 1/Ων0 with Ων0=0
        for l in 3:lmax
            y[i_Nl(l)] = c*k/((2*l+1)*aH(par,x0)) * y[i_Nl(l-1)]
        end
    else
        y[i_Nl(0):i_Nl(lmax)] .= 0.0
    end

    return y
end

#=
function perturbations_mode_tight(rec::Recombination, k; x1=-20.0, x2=0.0, kwargs...)
    function dy_dx!(x, y, dy)
        # pre-compute some common combined quantities
        ck_aH = c*k / aH(co,x)
        daH_aH = daH(co,x) / aH(co,x)
        R = 4*co.Ωγ0 / (3*co.Ωb0*a(x))

        # 1) un-pack variables from vector
        δc = y[i_δc]
        vc = y[i_vc]
        δb = y[i_δb]
        vb = y[i_vb]
        Φ  = y[i_Φ]
        N0 = y[i_Nl(0)]
        Nl = @view y[i_Nl(1):i_Nl(lmax)]
        Θ0 = y[i_Θl(0)]
        Θ1 = y[i_Θl(1)]
        Θ2 = (polarization ? -8/15 : -20/45) * ck_aH * Θ1/dτ(co,x)

        # 2) compute derivatives
        Ψ   = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θ2 + co.Ων0*Nl[2])
        dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0 + 4*co.Ων0/a(x)^2*N0)
        dδc = ck_aH*vc - 3*dΦ
        dδb = ck_aH*vb - 3*dΦ
        dvc = -vc - ck_aH*Ψ
        dΘ0 = -ck_aH*Θ1 - dΦ

        # calculate (dΘ1, dΘ2) (and dvb) with fixed-point iteration, as suggested by Callin
        d_aHdτ = daH(co,x)*dτ(co,x) + aH(co,x)*d2τ(co,x) # needed in fixed-point iteration
        function dvb_dΘ1_dΘ2_fixed_point(dvb_dΘ1_dΘ2)
            _, dΘ1, dΘ2 = dvb_dΘ1_dΘ2 # ignore dvb; it's in the tuple just so it is available after the fixed-point iteration
            d_3Θ1_plus_vb = (-((1-R)*dτ(co,x)+(1+R)*d2τ(co,x))*(3*Θ1+vb) - ck_aH*Ψ + (1-daH_aH)*ck_aH*(-Θ0+2*Θ2) + ck_aH*(-dΘ0+2*dΘ2)) / # notes approximate Θ2′=0 here (see discussion under Callin (34))
                            ((1+R)*dτ(co,x) + daH_aH - 1) # Callin's 3*Θ1′ + vb′ = Hans' q
            dvb = 1/(1+R) * (-vb - ck_aH*Ψ + R*(d_3Θ1_plus_vb+ck_aH*(-Θ0+2*Θ2) - ck_aH*Ψ))
            dΘ1 = (d_3Θ1_plus_vb - dvb) / 3
            dΘ2 = (polarization ? -8/15 : -20/45) * ck_aH * (dΘ1/dτ(co,x) - Θ1*d_aHdτ / (aH(co,x)*dτ(co,x)^2)) # anal Θ2′(x) # TODO: or from q?
            return (dvb, dΘ1, dΘ2)
        end
        dvb, dΘ1, _ = fixed_point_iterate(dvb_dΘ1_dΘ2_fixed_point, (NaN, 0.0, 0.0); tol=1e-10) # don't need dΘ2 below

        # re-pack variables into vector
        dy[i_δc]    = dδc
        dy[i_vc]    = dvc
        dy[i_δb]    = dδb
        dy[i_vb]    = dvb
        dy[i_Φ]     = dΦ
        dy[i_Θl(0)] = dΘ0
        dy[i_Θl(1)] = dΘ1

        # neutrinos
        if neutrinos
            dy[i_Nl(0)] = -ck_aH*Nl[1] - dΦ
            dy[i_Nl(1)] =  ck_aH/3 * (N0 - 2*Nl[2] + Ψ)
            for l in 2:lmax-1
                dy[i_Nl(l)] = ck_aH/(2*l+1) * (l*Nl[l-1] - (l+1)*Nl[l+1])
            end
            dy[i_Nl(lmax)] = ck_aH*Nl[lmax-1] - (lmax+1)/(aH(co,x)*η(co,x)) * Nl[lmax]
        else
            dy[i_Nl(0):i_Nl(lmax)] .= 0.0
        end

        return nothing
    end

    splines = Vector{Spline1D}(undef, i_max_full)
    y1 = perturbations_initial_conditions(co, x1, k)[1:i_max_tight] # cut away variables not in the tight system
    x, splines_tight = _spline_integral(dy_dx!, x1, x2, y1; name="perturbations tight (k=$(k*Mpc)/Mpc)", kwargs...)
    splines[1:i_max_tight] .= splines_tight

    # extend tight splines to full system
    splines[i_Θl(2)] = spline(x, @. (polarization ? -8/15 : -20/45) * c*k / (aH(co,x)*dτ(co,x)) * splines[i_Θl(1)](x))
    for l in 3:lmax
        splines[i_Θl(l)] = spline(x, @. -l/(2*l+1) * c*k/(aH(co,x)*dτ(co,x)) * splines[i_Θl(l-1)](x))
    end

    if polarization
        splines[i_ΘPl(0)] = spline(x, @. 5/4 * splines[i_Θl(2)](x))
        splines[i_ΘPl(1)] = spline(x, @. -c*k/(4*aH(co,x)*dτ(co,x)) * splines[i_Θl(2)](x))
        splines[i_ΘPl(2)] = spline(x, @. 1/4 * splines[i_Θl(2)](x))
        for l in 3:lmax
            splines[i_ΘPl(l)] = spline(x, @. -l/(2*l+1) * c*k/(aH(co,x)*dτ(co,x)) * splines[i_ΘPl(l-1)](x))
        end
    else
        for l in 0:lmax
            splines[i_ΘPl(l)] = spline(x, 0 .* x)
        end
    end

    return x, splines
end
=#

# putting splines for η and τ behind a function barrier removes allocations!
function perturbations_mode_full(par::Parameters, η::ODESolution, τ::ODESolution, k; x1=-20.0, x2=0.0)
    function dy_dx!(dy, y, _, x)
        # 1) un-pack variables from vector
        δc = y[i_δc]
        vc = y[i_vc]
        δb = y[i_δb]
        vb = y[i_vb]
        Φ  = y[i_Φ]
        Nl = OffsetVector(view(y, i_Nl(0):i_Nl(lmax)), 0:lmax)
        Θl = OffsetVector(view(y, i_Θl(0):i_Θl(lmax)), 0:lmax)
        ΘPl = OffsetVector(view(y, i_ΘPl(0):i_ΘPl(lmax)), 0:lmax)

        # 2) pre-compute some commonly combined auxilliary quantities
        ck_aH = c*k / aH(par,x)
        R = 4*par.Ωγ0 / (3*par.Ωb0*a(x))
        τ′x = ForwardDiff.derivative(τ, x)
        ηx = η(x)

        # 3) compute derivatives
        Ψ   = -Φ - 12 * (par.H0 / (c*k*a(x)))^2 * (par.Ωγ0*Θl[2]+par.Ων0*Nl[2])
        Π   = Θl[2] + ΘPl[0] + ΘPl[2]
        dΦ  = Ψ - ck_aH^2/3*Φ + (par.H0/aH(par,x))^2/2 * (par.Ωc0/a(x)*δc + par.Ωb0/a(x)*δb + 4*par.Ωγ0/a(x)^2*Θl[0] + 4*par.Ων0/a(x)^2*Nl[0])
        dδc = ck_aH*vc - 3*dΦ
        dδb = ck_aH*vb - 3*dΦ
        dvc = -vc - ck_aH*Ψ
        dvb = -vb - ck_aH*Ψ + τ′x*R*(3*Θl[1]+vb)
        dΘl0   = -ck_aH*Θl[1] - dΦ
        dΘl1   =  ck_aH/3 * (Θl[0]-2*Θl[2]+Ψ) + τ′x*(Θl[1]+vb/3)
        dΘlmax = ck_aH*Θl[lmax-1] - (lmax+1)/(aH(par,x)*ηx)*Θl[lmax] + τ′x*Θl[lmax] # 2nd term: their η is my c*η
        dΘPl0   = -ck_aH*ΘPl[1] + τ′x*(ΘPl[0]-Π/2)
        dΘPlmax =  ck_aH*ΘPl[lmax-1] - (lmax+1)/(aH(par,x)*ηx)*ΘPl[lmax] + τ′x*ΘPl[lmax]

        # 4) re-pack variables into vector
        dy[i_δc]    = dδc
        dy[i_vc]    = dvc
        dy[i_δb]    = dδb
        dy[i_vb]    = dvb
        dy[i_Φ]     = dΦ
        dy[i_Θl(0)] = dΘl0
        dy[i_Θl(1)] = dΘl1
        for l in 2:lmax-1
            dy[i_Θl(l)] = ck_aH/(2*l+1) * (l*Θl[l-1] - (l+1)*Θl[l+1]) + τ′x*(Θl[l]-Π/10*δ(l,2))
        end
        dy[i_Θl(lmax)] = dΘlmax

        if polarization
            dy[i_ΘPl(0)] = dΘPl0
            for l in 1:lmax-1
                dy[i_ΘPl(l)] = ck_aH/(2*l+1) * (l*ΘPl[l-1] - (l+1)*ΘPl[l+1]) + τ′x*(ΘPl[l]-Π/10*δ(l,2))
            end
            dy[i_ΘPl(lmax)] = dΘPlmax
        else
            dy[i_ΘPl(0):i_ΘPl(lmax)] .= 0.0
        end

        if neutrinos
            dy[i_Nl(0)] = -ck_aH*Nl[1] - dΦ
            dy[i_Nl(1)] =  ck_aH/3 * (Nl[0] - 2*Nl[2] + Ψ)
            for l in 2:lmax-1
                dy[i_Nl(l)] = ck_aH/(2*l+1) * (l*Nl[l-1] - (l+1)*Nl[l+1])
            end
            dy[i_Nl(lmax)] = ck_aH*Nl[lmax-1] - (lmax+1)/(aH(par,x)*ηx) * Nl[lmax]
        else
            dy[i_Nl(0):i_Nl(lmax)] .= 0.0
        end

        return nothing # dy is in-place
    end

    t1 = now()
    y1 = perturbations_initial_conditions(par, η, τ, x1, k)
    prob = ODEProblem{true}(dy_dx!, y1, (x1, x2))
    sol = solve(prob, KenCarp4(autodiff=false); abstol=1e-8, reltol=1e-8)
    t2 = now()
    println("Integrated perturbation mode k=$(k*Mpc)/Mpc with $(length(sol.t)) points in $(t2-t1)")
    return sol
end

function perturbations_mode(rec::Recombination, k; tight::Bool=false, kwargs...)
    #=
    bg = rec.bg
    par = bg.par

    splines = Vector{Spline1D}(undef, i_max_ext)
    if tight
        # 1) use tight coupling approximation at early times to avoid stiff equations,
        # 2) then integrate the full equations from when the approximation breaks down,
        # 3) merge 1+2
        x12 = time_tight_coupling(co, k)
        x1, spl1s = perturbations_mode_tight(co, k; x2=x12, solver=Tsit5(), kwargs...)
        y12 = [spl1(x12) for spl1 in spl1s] # give final tight values as ICs for full system
        x2, spl2s = perturbations_mode_full(co, k; x1=x12, y1=y12, solver=Tsit5(), kwargs...)
        x, spls = splinejoin(x1, x2, spl1s, spl2s)
        splines[1:i_max_full] .= spls
    else
    =#
        # only integrate the full (stiff) equations using an appropriate solver
        # discussion of stiff solvers / tight coupling etc. in context of Boltzmann solvers / Julia / DifferentialEquations:
        # https://discourse.julialang.org/t/is-autodifferentiation-possible-in-this-situation/54807
        return perturbations_mode_full(rec.bg.par, rec.bg.η, rec.τ, k) # KenCarp4(autodiff=false) and radau() work well!  # TODO: autodiff here
        #x, splines_full = perturbations_mode_full(rec, k; solver=KenCarp4(autodiff=false), abstol=1e-10, reltol=1e-10, verbose=true, kwargs...) # KenCarp4(autodiff=false) and radau() work well!
        #splines[1:i_max_full] .= splines_full
    #end

    # extend with variables given in terms of the integrated ones (like Ψ and S)
    #=
    splines[i_Ψ] = spline(x, @. -splines[i_Φ](x) - 12*par.H0^2/(c*k*a(x))^2 * (par.Ωγ0*splines[i_Θl(2)](x) + par.Ων0*splines[i_Nl(2)](x)))

    Π            = spline(x, @.  splines[i_Θl(2)](x) + splines[i_ΘPl(0)](x) + splines[i_ΘPl(2)](x))
    aH_g_vb      = spline(x, @. aH(par,x) * g(rec,x) * splines[i_vb](x))
    aH_g_Π       = spline(x, @. aH(par,x) * g(rec,x) * Π(x))
    aH_d_aH_g_Π  = spline(x, aH.(par,x) .* derivative(aH_g_Π, x))
    splines[i_S] = spline(x, g.(rec,x) .* (splines[i_Θl(0)](x) .+ splines[i_Ψ](x) .+ Π(x)/4) .+
                             exp.(-τ.(rec,x)) .* (derivative(splines[i_Ψ], x) .- derivative(splines[i_Φ], x)) .-
                             1/(c*k) * derivative(aH_g_vb, x) .+
                             3/(4*c^2*k^2) * derivative(aH_d_aH_g_Π, x)
                   )
    return x, splines
    =#
end

#=
function perturbations_splines(co::ΛCDM; tight::Bool=false, ks=nothing)
    if length(co.perturbation_splines2D) == 0 || !isnothing(ks)
        t1 = now()

        if isnothing(ks)
            kmin =    1 / (c*η(co,0))
            kmax = 3000 / (c*η(co,0))
            ks = 10 .^ (log10(kmin) .+ (log10(kmax)-log10(kmin)) * range(0, 1; length=200))
        end
        
        # (x,k) spline requires 1D arrays for x and k (and does not accept a fully irregular 2D grid)
        # use the x-values from the mode that has the most (not necessarily k=kmax)
        spliness = Vector{Vector{Spline1D}}(undef, length(ks))
        xs = nothing
        @threads for i_k in 1:length(ks) # call julia with --threads=8 to get a decent speed-up
            k = ks[i_k]
            x, spliness[i_k] = perturbations_mode(co, k; tight=tight)
            if isnothing(xs) || length(x) > length(xs)
                xs = x
            end
        end

        perturbs = Array{Float64, 3}(undef, i_max_ext, length(xs), length(ks)) # indexed as [i_quantity, i_x, i_k]
        for i_k in 1:length(ks)
            for i_qty in 1:i_max_ext
                perturbs[i_qty, :, i_k] .= spliness[i_k][i_qty](xs)
            end
        end

        co.perturbation_splines2D = [] # reset if not empty
        for i_qty in 1:i_max_ext
            qty = @view perturbs[i_qty, :, :]
            push!(co.perturbation_splines2D, Spline2D(xs, ks, qty)) # spline (x, k)
        end

        t2 = now()
        println("Splined perturbations(x,k) on ($(length(xs)),$(length(ks))) grid in $(t2-t1) using $(sizeof(perturbs)/1e6) MB and $(nthreads()) parallel threads")
    end
    return co.perturbation_splines2D
end
=#

#=
function perturbations_mode_splines(co::ΛCDM, k; kwargs...)
    i = searchsortedfirst(co.perturbation_splines1D, k; by = tuple -> tuple[1]) # index of existing k-mode
    if i > length(co.perturbation_splines1D) || co.perturbation_splines1D[i][1] != k # not computed before
        _, splines = perturbations_mode(co, k; kwargs...)
        insert!(co.perturbation_splines1D, i, (k, splines))
    end
    return co.perturbation_splines1D[i][2]
end
=#

# raw quantities (from integration)
Φ(perts::PerturbationMode, x) = perts.y(x)[i_Φ]
δc(perts::PerturbationMode, x) = perts.y(x)[i_δc]
δb(perts::PerturbationMode, x) = perts.y(x)[i_δb]
vc(perts::PerturbationMode, x) = perts.y(x)[i_vc]
vb(perts::PerturbationMode, x) = perts.y(x)[i_vb]
Θl(perts::PerturbationMode, x, l) = perts.y(x)[i_Θl(l)]
Nl(perts::PerturbationMode, x, l) = perts.y(x)[i_Nl(l)]
ΘPl(perts::PerturbationMode, x, l) = perts.y(x)[i_ΘPl(l)]

# "composite" quantities (form raw quantities)
function Ψ(perts::PerturbationMode, x; deriv=0)
    par = perts.rec.bg.par
    k = perts.k
    return -Φ(perts, x) - 12*par.H0^2/(c*k*a(x))^2 * (par.Ωγ0*Θl(perts, x, 2) + par.Ων0*Nl(perts, x, 2))
end

Π(perts::PerturbationMode, x) = Θl(perts, x, 2) + ΘPl(perts, x, 0) + ΘPl(perts, x, 2)

# TODO: rename S -> ST?
S_SW(perts::PerturbationMode, x) = g(perts.rec,x) * (Θl(perts, x, 0) + Ψ(perts, x) + Π(perts,x)/4)
S_ISW(perts::PerturbationMode, x) = exp(-τ(perts.rec,x)) * ForwardDiff.derivative(x -> Ψ(perts,x) - Φ(perts,x), x)
S_Doppler(perts::PerturbationMode, x) = -1/(c*perts.k) * ForwardDiff.derivative(x -> aH(perts.rec.bg.par,x) * g(perts.rec,x) * vb(perts, x), x)
S_polarization(perts::PerturbationMode, x) = 3/(4*c^2*perts.k^2) * ForwardDiff.derivative(x -> aH(perts.rec.bg.par,x) * ForwardDiff.derivative(x -> aH(perts.rec.bg.par,x) * g(perts.rec,x) * Π(perts,x), x), x)
S(perts::PerturbationMode, x) = S_SW(perts, x) + S_ISW(perts, x) + S_Doppler(perts, x) + S_polarization(perts, x)
SE(perts::PerturbationMode, x) = x == 0.0 ? 0.0 : 3 * g(perts.rec,x) * Π(perts,x) / (2*c*perts.k*(η(perts.rec.bg,0.0)-η(perts.rec.bg,x)))^2 # force SE(x=0) = 0 (otherwise gives Inf)
