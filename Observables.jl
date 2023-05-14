P_primordial(par::Parameters, k::Real) = 2*π^2/k^3 * par.As * (k/par.k_pivot)^(par.ns-1)
Δ(perturb::PerturbationMode, x::Real, k::Real) = c^2*k^2*Φ(perturb,x) / (3/2*perturb.rec.bg.par.Ωm0/a(x)*perturb.rec.bg.par.H0^2)
P(perturb::PerturbationMode, x::Real, k::Real) = abs(Δ(perturb,x,k))^2 * P_primordial(perturb.rec.bg.par, k)

# TODO:
# From Bolt.jl creator: Levin method for oscillatory Bessel integrals
# https://discourse.julialang.org/t/rfc-ann-oscillatoryintegralsode-jl-levin-method-ordinarydiffeq/55601

function integrate_trapz(f, x1, x2; step=nothing, length=nothing)
    xs = range(x1, x2; step=step, length=length) # either step or length should be specified
    return trapz(xs, f.(xs))
end

function integrate_adaptive(f, x1, x2; atol=0, rtol=1e-3, order=8) # TODO: rtol=1e-3 and order=8 gives good results!
    return quadgk(f, x1, x2; atol=atol, rtol=rtol, order=order)[1] # discard error
end

#integratex_default = (f, x1, x2) -> integrate_trapz(f, x1, x2; length=1000)
#integratek_default = (f, k1, k2) -> integrate_trapz(f, k1, k2; step=1000)
function Cl(l, S, η, par; )
    η0 = η(0.0)

    integratex = integrate_adaptive
    integratek = integrate_adaptive

    dΘ0_dx(x, k) :: Float64 = S(x, k) * jl(l, c*k*(η0-η(x)))
    Θ0(k) = integratex(x -> dΘ0_dx(x, k), -20.0, 0.0)

    dCl_dk(k) :: Float64 = 2/π * P_primordial(par, k) * (k * Θ0(k))^2
    Cl = integratek(dCl_dk, 1/(c*η0), 4000/(c*η0))

    return Cl
end

#=
function Cl_trapz(l, Sspl, η, par, Ny_per_osc=20, Nk_per_osc=10)
    η0 = η(0.0)

    # TODO: allocate x-array ONCE

    Ny = 1000
    dΘ0_dk_grid = Vector{Float64}(undef, Ny)
    function Θ0(k)
        # TODO: set ymin, ymax instead?
        #xmin = -20.0
        #xmax = 0.0
        #y(x) = c*k*(η0-η(x))
        #ymin = y(xmax)
        #ymax = y(xmin)
        #Δy = 2*π/Ny_per_osc
        #Ny = Int(round((ymax-ymin) / Δy))
        #Ny = max(Ny, 1000) # TODO: do the maximum with what is needed to resolve source function
        #Ny = 1000
        xs = range(xmin, xmax, length=Ny)
        ys = y.(xs)
        dΘ0_dk_grid .= Sspl.(xs, k) .* jl.(l, c*k * (η0 .- η.(xs))) # (x,k) function
        return trapz(xs, dΘ0_dk_grid) # integrate over x
    end

    # TODO: spline jl. all allocations are from it!

    # TODO: @code_warntype gives η(xs)::Any ???!?!? use η.(xs) instead?
    # TODO: jl allocates a lot. reduce somehow?
    #dΘ0_dk_grid = evalgrid(Sspl, xs, ks) .* jl.(l, c * ks' .* (η(0.0) .- η.(xs))) # (x,k) function

    kmin =    1 / (c*η(0))
    kmax = 4000 / (c*η(0))
    Δk = 2*π / (c*η0*Nk_per_osc)
    ks = range(kmin, kmax, step=Δk)

    println("$(length(ks)) k-values")

    dCl_dk_grid = 2/π * P_primordial.(par, ks) .* (ks .* Θ0.(ks)) .^ 2
    return trapz(ks, dCl_dk_grid) # integrate over k
end
=#

# TODO: for some reason it fails for l=8 and rtol=1e-1
#=
function Cl_cubature(l, Sspl, η, par; rtol=1e-3, kwargs...)
    η0 = η(0.0)
    S(x,k) = Sspl(x, log10(k)) # TODO: log or not?
    function integrand(vec)
        k, x1, x2 = vec[1], vec[2], vec[3]
        return S(x1, k) * jl(l, c*k*(η0 - η(x1))) * # dΘ0_dx (#1)
               S(x2, k) * jl(l, c*k*(η0 - η(x2))) * # dΘ0_dx (#2)
               2/π * k^2 * P_primordial(par, k)   # dCl_dk
    end
    kmin = 1 / (c*η0)
    kmax = 4000 / (c*η0)
    xmin = -20.0
    xmax = 0.0
    return hcubature(integrand, (kmin, xmin, xmin), (kmax, xmax, xmax); rtol=rtol, kwargs...)[1]
end
=#

function spline_S(rec)
    kmin =    1 / (c*η(rec.bg,0))
    kmax = 4000 / (c*η(rec.bg,0))
    logks = range(log10(kmin), log10(kmax), length=250) # TODO: 250
    ks = 10 .^ logks

    perturbs = [PerturbationMode(rec, k) for k in ks]
    #index_with_most_points = argmax(length(perturb.qty_splines.t) for perturb in perturbs)
    #xs = perturbs[index_with_most_points].y.t # TODO: extend?

    # make uniform grids
    xs = range(-20, 0, length=2000) # TODO: looks like this is enough?

    # TODO: evaluate quicker
    Sdata = Float64[S(perturb, x) for x in xs, perturb in perturbs]

    println("Splining S(x, log10(k)) on uniform $(size(Sdata))-grid")
    S_spline_x_logk = spline((xs, logks), Sdata) # (x, log10(k)) spline - must convert to S(x,k) spline upon calling it!
    S_spline_x_k(x, k) = S_spline_x_logk(x, log10(k))
    return S_spline_x_k
end

function Cls(rec::Recombination, ls::Vector)
    bg = rec.bg
    S = spline_S(rec) # (x,log(k))

    # TODO: paralellize with threads?
    Clarr = Vector{Float64}(undef, length(ls))
    Threads.@threads for i in 1:length(ls)
        l = ls[i]
        @time(time = @elapsed(Clarr[i] = Cl(l, S, bg.η, bg.par)))
        println("C(l=$l) = $(Clarr[i]) ($time seconds)")
    end

    return Clarr
end
