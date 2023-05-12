P_primordial(par::Parameters, k::Real) = 2*π^2/k^3 * par.As * (k/par.k_pivot)^(par.ns-1)
Δ(perturb::Perturbations, x::Real, k::Real) = c^2*k^2*Φ(perturb,x) / (3/2*perturb.rec.bg.par.Ωm0/a(x)*perturb.rec.bg.par.H0^2)
P(perturb::Perturbations, x::Real, k::Real) = abs(Δ(perturb,x,k))^2 * P_primordial(perturb.rec.bg.par, k)

# From Bolt.jl creator: Levin method for oscillatory Bessel integrals
# https://discourse.julialang.org/t/rfc-ann-oscillatoryintegralsode-jl-levin-method-ordinarydiffeq/55601

# seems to work, but allocates a lot!
# TODO: instead of nested quadgk, use dedicated 2D integral algorithm? Cubature.jl?

# adaptive trapezoid method
function integrate(f, x1, x2, Δx)
    x = (x1 + x2) / 2
    if x2 - x1 < Δx(x)
        f1 = f(x1)
        f2 = f(x2)
        return 1/2 * (x2 - x1) * (f(x1) + f(x2)) # trapezoid area
    else
        return integrate(f, x1, x, Δx) + integrate(f, x, x2, Δx)
    end
end

function Cl_fancy(l::Integer, Sspl, η, par::Parameters)
    η0 = η(0.0)
    S(x, k) = Sspl(x, log10(k)) # TODO: needed?
    dΘ0_dx(x,k) :: Float64 = S(x, k) * jl(l, c*k*(η0 - η(x)))
    Θ0(k) = integrate(x -> dΘ0_dx(x,k), -20.0, 0.0, x -> 2*π*aH(par,x) / (10*c*k))
    dCl_dk(k) :: Float64 = 2/π * k^2 * P_primordial(par, k) * Θ0(k)^2
    return integrate(dCl_dk, 1/(c*η0), 3000/(c*η0), k -> 2*π / (10*c*η0))
end

function Cl_quadgk(l::Integer, Sspl, η, par::Parameters; rtol=1e-3, order=8, kwargs...)
    η0 = η(0.0)
    dΘ0_dx(x, k) ::Float64 = Sspl(x, k) * jl(l, c*k*(η0 - η(x)))
    Θ0(k) = quadgk(x -> dΘ0_dx(x,k), -20.0, 0.0; rtol=rtol, order=order, kwargs...)[1]
    dCl_dk(k) ::Float64 = 2/π * k^2 * P_primordial(par, k) * Θ0(k)^2
    return quadgk(dCl_dk, 1/(c*η0), 3000/(c*η0); rtol=rtol, order=order, kwargs...)[1] # TODO: dynamic kmin
end

# TODO: type-stable, little-allocating trapz version
function Cl_trapz(l, Sspl, η, par)
    kmin =    1 / (c*η(0))
    kmax = 3000 / (c*η(0))
    Δk = 2*π / (c*η(0)*20)
    ks = range(kmin, kmax, step=Δk)
    xs = range(-20.0, 0.0, length=1000) # TODO: take from perturbation solution? TODO: see Callin eq (51) for spacing

    # TODO: spline jl. all allocations are from it!

    # TODO: @code_warntype gives η(xs)::Any ???!?!? use η.(xs) instead?
    # TODO: jl allocates a lot. reduce somehow?
    #dΘ0_dk_grid = evalgrid(Sspl, xs, ks) .* jl.(l, c * ks' .* (η(0.0) .- η.(xs))) # (x,k) function
    dΘ0_dk_grid = Sspl.(xs, log10.(ks)') .* jl.(l, c * ks' .* (η(0.0) .- η.(xs))) # (x,k) function
    Θ0 = trapz(xs, dΘ0_dk_grid, Val(1)) # integrate over x
    dCl_dk_grid = 2/π * P_primordial.(par, ks) .* (ks .* Θ0) .^ 2
    return trapz(ks, dCl_dk_grid) # integrate over k
end

# TODO: for some reason it fails for l=8 and rtol=1e-1
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
    kmax = 3000 / (c*η0)
    xmin = -20.0
    xmax = 0.0
    return hcubature(integrand, (kmin, xmin, xmin), (kmax, xmax, xmax); rtol=rtol, kwargs...)[1]
end

function S_spline(rec)
    kmin =    1 / (c*η(rec.bg,0))
    kmax = 3000 / (c*η(rec.bg,0))
    logks = range(log10(kmin), log10(kmax), length=250)
    ks = 10 .^ logks

    perturbs = [Perturbations(rec, k) for k in ks]
    xs = perturbs[argmax(length(perturb.qty_splines.t) for perturb in perturbs)].qty_splines.t # TODO: extend?

    # make uniform grids
    xs = range(-20, 0, length=5000) # TODO: looks like this is enough?

    # TODO: evaluate quicker
    Sdata = Float64[S(perturb, x) for x in xs, perturb in perturbs]

    println("Splining S(x, log10(k)) on uniform $(size(Sdata))-grid")
    S_spline_x_logk = scale(interpolate(Sdata, BSpline(Cubic(Line(OnGrid())))), (xs, logks)) # (x, log10(k)) spline - must convert to S(x,k) spline upon calling it!

    S_spline_x_k(x, k) = S_spline_x_logk(x, log10(k))
    return S_spline_x_k
end

function Cls(rec::Recombination, ls::Vector)
    bg = rec.bg
    Sspl = S_spline(rec) # (x,log(k))

    # TODO: paralellize with threads?
    Clarr = Vector{Float64}(undef, length(ls))
    for (i, l) in enumerate(ls)
        time = @elapsed(Clarr[i] = Cl_quadgk(l, Sspl, bg.η_spline, bg.par))
        println("C(l=$l) = $(Clarr[i]) ($time seconds)")
    end

    return Clarr
end
