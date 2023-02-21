module Cosmology

export a, z
export ΛCDM
export H, aH, daH, d2aH
export t, η, Ωγ, Ων, Ωb, Ωc, Ωk, ΩΛ, Ωr, Ωm, Ω
export r_m_equality, m_Λ_equality
export dL

include("Constants.jl")

using DifferentialEquations
using Interpolations
using Roots
using .Constants

struct ΛCDM
    h::Float64 # dimensionless Hubble parameter h = H0 / (100 TODO units)
    H0::Float64
    Ωb0::Float64 # baryons
    Ωc0::Float64 # cold dark matter
    Ωm0::Float64 # matter = baryons + cold dark matter
    Ωk0::Float64 # curvature
    Ωγ0::Float64 # photons
    Ων0::Float64 # (massless) neutrinos
    Ωr0::Float64 # radiation = photons + neutrinos
    ΩΛ::Float64 # dark energy (as cosmological constant)
    Tγ0::Float64 # photon (CMB) temperature
    Neff::Float64 # TODO: ?

    x_spline::Vector{Float64}
    η_spline::Interpolations.Extrapolation # TODO: what kind of spline?
    t_spline::Interpolations.Extrapolation

    function ΛCDM(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046)
        H0 = h * 100*km/Mpc
        Ωm0 = Ωb0+Ωc0
        Ωγ0 = 2 * π^2/30 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ  = 1.0 - (Ωr0 + Ωm0 + Ωk0)

        Ω0 = Ωr0 + Ωm0 + Ωk0 + ΩΛ

        x1, x2 = -20, +20 # integration and spline range
        Hx(x::Real) = H(h, Ωr0, Ωm0, Ωk0, ΩΛ, x)

        dη_dx(η, p, x) = c / (a(x) * Hx(x)) # TODO: integrate in dimensionless units closer to 1
        η1 = c / (a(x1) * Hx(x1))
        condition(η, x, integrator) = Ωpoly(Ωr0, Ωm0, Ωk0, ΩΛ, x) - 1e-4
        affect!(integrator) = terminate!(integrator)
        big_rip_terminator = ContinuousCallback(condition, affect!)
        sol = solve(ODEProblem(dη_dx, η1, (x1, x2)), Tsit5(); reltol=1e-10, callback = big_rip_terminator) # automatically choose method
        x, η = sol.t, sol.u
        if x[end] == x[end-1]
            # remove last duplicate point TODO: a better way?
            x = x[1:end-1]
            η = η[1:end-1]
        end
        η_spline = linear_interpolation(x, η) # spline

        dt_dx(t, p, x) = 1 / Hx(x)
        t1 = 1 / (2*Hx(x1))
        sol = solve(ODEProblem(dt_dx, t1, (x1, x2)), Tsit5(); reltol=1e-10, callback = big_rip_terminator)
        x, t = sol.t, sol.u
        if x[end] == x[end-1]
            # remove last duplicate point TODO: a better way?
            x = x[1:end-1]
            t = t[1:end-1]
        end
        t_spline = linear_interpolation(x, t) # spline

        # TODO: don't pass x: it is different from η to t
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ, Tγ0, Neff, x, η_spline, t_spline)
    end
end

# TODO: how to not vectorize over the cosmology without passing it as a keyword argument? https://docs.julialang.org/en/v1/manual/functions/#man-vectorized
# TODO: relevant discussion: https://discourse.julialang.org/t/recent-broadcast-changes-iterate-by-default-scalar-struct-and/11178/26
Base.broadcastable(co::ΛCDM) = Ref(co) # never broadcast cosmology (in vectorized calls)

# internal time variable: natural logarithm of scale factor
x(a::Real) = log(a)
a(x::Real) = exp(x)
z(x::Real) = 1/a(x) - 1

# Friedmann equation
Ωpoly(Ωr0::Real, Ωm0::Real, Ωk0::Real, ΩΛ::Real, x::Real; d::Integer=0) = (-4)^d * Ωr0/a(x)^4 + (-3)^d * Ωm0/a(x)^3 + (-2)^d * Ωk0/a(x)^2 + 0^d * ΩΛ # evolution of densities relative to *today's* critical density # TODO: can do derivative here (0^0 = 1 in Julia)! TODO: call it Ωpoly?
Ωpoly(co::ΛCDM, x::Real; d::Integer=0) = Ωpoly(co.Ωr0, co.Ωm0, co.Ωk0, co.ΩΛ, x; d)
#dΩevo(Ωr0::Real, Ωm0::Real, Ωk0::Real, ΩΛ0::Real, x::Real; n::Integer=1) = (-4)^n * Ωr0/a(x)^4 + (-3)^n * Ωm0/a(x)^3 + (-2)^n * Ωk0/a(x)^2 # n-th (n >= 1) derivative of Ωev0 wrt. x
#dΩevo(co::ΛCDM, x::Real; n=Integer::1) = dΩevo(co.Ωr0, co.Ωm0, co.Ωk0, co.ΩΛ0, x; n)
H(h::Real, Ωr0::Real, Ωm0::Real, Ωk0::Real, ΩΛ::Real, x::Real) = h * 100 * km/Mpc * √(max(Ωpoly(Ωr0, Ωm0, Ωk0, ΩΛ, x), 0.0)) # in case we have a slight negative value due to floating point arithmetic, add 0im and take real part
H(co::ΛCDM, x::Real) = H(co.h, co.Ωγ0+co.Ων0, co.Ωb0+co.Ωc0, co.Ωk0, co.ΩΛ, x) # "normal" Hubble parameter

# conformal Hubble parameter (𝓗 = a*H) + 1st derivative + 2nd derivative
aH(co::ΛCDM, x::Real) = a(x) * H(co, x) # conformal Hubble parameter
daH(co::ΛCDM, x::Real) = aH(co, x) * (1 + 1/2 * Ωpoly(co, x; d=1) / Ωpoly(co, x))
#d2aH(co::ΛCDM, x::Real) = daH(co, x)^2 / ah(co, x) + 1/2 * aH(co, x) * (Ωpoly(co, x; d=2) / Ωpoly(co, x) + (Ωpoly(co, x; d=1) / Ωpoly(co, x))^2)
d2aH(co::ΛCDM, x::Real) = aH(co, x) * (1 + Ωpoly(co, x; d=1) / Ωpoly(co, x) + 1/2 * Ωpoly(co, x; d=2) / Ωpoly(co, x) - 1/4 * (Ωpoly(co, x; d=1) / Ωpoly(co, x))^2)

# conformal time + cosmic time
η(co::ΛCDM, x::Real) = co.η_spline(x)
t(co::ΛCDM, x::Real) = co.t_spline(x)

# density parameters (relative to critical density *at the time*)
Ωγ(co::ΛCDM, x::Real) = co.Ωγ0 / (a(x)^4 * H(co, x)^2 / H(co, 0)^2)
Ων(co::ΛCDM, x::Real) = co.Ων0 / (a(x)^4 * H(co, x)^2 / H(co, 0)^2)
Ωb(co::ΛCDM, x::Real) = co.Ωb0 / (a(x)^3 * H(co, x)^2 / H(co, 0)^2)
Ωc(co::ΛCDM, x::Real) = co.Ωc0 / (a(x)^3 * H(co, x)^2 / H(co, 0)^2)
Ωk(co::ΛCDM, x::Real) = co.Ωk0 / (a(x)^2 * H(co, x)^2 / H(co, 0)^2)
ΩΛ(co::ΛCDM, x::Real) = co.ΩΛ  / (H(co, x)^2 / H(co, 0)^2)
Ωr(co::ΛCDM, x::Real) = Ωγ(co, x) + Ων(co, x)
Ωm(co::ΛCDM, x::Real) = Ωb(co, x) + Ωc(co, x)
Ω(co::ΛCDM, x::Real) = Ωr(co, x) + Ωm(co, x) + Ωk(co, x) + ΩΛ(co, x)

# equalities
r_m_equality(co::ΛCDM) = find_zero(x -> Ωr(co, x) - Ωm(co, x), (-20, +20)) # TODO: save xmin, xmax in ΛCDM
m_Λ_equality(co::ΛCDM) = find_zero(x -> Ωm(co, x) - ΩΛ(co, x), (-20, +20)) # TODO: save xmin, xmax in ΛCDM

# distances
function dL(co::ΛCDM, x::Real)
    χ = η(co, 0) - η(co, x)

    # TODO: change back
    # Ωk0 = 0: r = χ
    # Ωk0 > 0: r = χ *  sinc(√(Ωk0)*H0*χ/c)
    # Ωk0 < 0: r = χ * sinhc(√(Ωk0)*H0*χ/c)
    # ... are all captured by
    r = χ * real(sinc(√(complex(co.Ωk0)) * co.H0 * χ / c))
    #=
    if co.Ωk0 < 0
        r = χ *  sin(√(-co.Ωk0)*co.H0*χ/c) / (√(-co.Ωk0)*co.H0*χ/c)
    elseif co.Ωk0 > 0
        r = χ * sinh(√(+co.Ωk0)*co.H0*χ/c) / (√(+co.Ωk0)*co.H0*χ/c)
    else
        r = χ
    end
    =#


    return r / a(x)
end

end
