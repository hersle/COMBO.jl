module Cosmology

export a, z
export ΛCDM
export H, aH, daH, d2aH
export t, η, Ωγ, Ων, Ωb, Ωc, Ωk, ΩΛ, Ωr, Ωm, Ω
export r_m_equality, m_Λ_equality, acceleration_onset
export dL

include("Constants.jl")

using DifferentialEquations
using Dierckx # spline interpolation
using Roots
using .Constants

mutable struct ΛCDM
    const h::Float64 # dimensionless Hubble parameter h = H0 / (100 TODO units)
    const H0::Float64
    const Ωb0::Float64 # baryons
    const Ωc0::Float64 # cold dark matter
    const Ωm0::Float64 # matter = baryons + cold dark matter
    const Ωk0::Float64 # curvature
    const Ωγ0::Float64 # photons
    const Ων0::Float64 # (massless) neutrinos
    const Ωr0::Float64 # radiation = photons + neutrinos
    const ΩΛ::Float64 # dark energy (as cosmological constant)
    const Tγ0::Float64 # photon (CMB) temperature
    const Neff::Float64 # TODO: ?

    η_spline::Union{Nothing, Dierckx.Spline1D} # TODO: what kind of spline?
    t_spline::Union{Nothing, Dierckx.Spline1D}

    function ΛCDM(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046)
        H0 = h * 100*km/Mpc
        Ωm0 = Ωb0 + Ωc0
        Ωγ0 = 2 * π^2/30 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ  = 1.0 - (Ωr0 + Ωm0 + Ωk0)
        Ω0 = Ωr0 + Ωm0 + Ωk0 + ΩΛ # ≈ 1
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ, Tγ0, Neff, nothing, nothing)
    end
end

Base.broadcastable(co::ΛCDM) = Ref(co) # never broadcast cosmology (in vectorized calls)

# internal time variable: natural logarithm of scale factor
x(a::Real) = log(a)
a(x::Real) = exp(x)
z(x::Real) = 1/a(x) - 1

# Friedmann equation
# TODO: rename E
Ωpoly(Ωr0::Real, Ωm0::Real, Ωk0::Real, ΩΛ::Real, x::Real; d::Integer=0) = (-4)^d * Ωr0/a(x)^4 + (-3)^d * Ωm0/a(x)^3 + (-2)^d * Ωk0/a(x)^2 + 0^d * ΩΛ # evolution of densities relative to *today's* critical density # TODO: can do derivative here (0^0 = 1 in Julia)! TODO: call it Ωpoly?
Ωpoly(co::ΛCDM, x::Real; d::Integer=0) = Ωpoly(co.Ωr0, co.Ωm0, co.Ωk0, co.ΩΛ, x; d)
H(h::Real, Ωr0::Real, Ωm0::Real, Ωk0::Real, ΩΛ::Real, x::Real) = h * 100 * km/Mpc * √(Ωpoly(Ωr0, Ωm0, Ωk0, ΩΛ, x))
H(co::ΛCDM, x::Real) = H(co.h, co.Ωγ0+co.Ων0, co.Ωb0+co.Ωc0, co.Ωk0, co.ΩΛ, x) # "normal" Hubble parameter

# conformal Hubble parameter (𝓗 = a*H) + 1st derivative + 2nd derivative (from analytical considerations)
  aH(co::ΛCDM, x::Real) = a(x) * H(co, x)
 daH(co::ΛCDM, x::Real) = aH(co, x) * (1 + 1/2 * Ωpoly(co, x; d=1) / Ωpoly(co, x))
d2aH(co::ΛCDM, x::Real) = aH(co, x) * (1 + Ωpoly(co, x; d=1) / Ωpoly(co, x) + 1/2 * Ωpoly(co, x; d=2) / Ωpoly(co, x) - 1/4 * (Ωpoly(co, x; d=1) / Ωpoly(co, x))^2)

function _spline_dy_dx(co::ΛCDM, dy_dx::Function, x1::Float64, x2::Float64, y1::Float64)
    sol = solve(ODEProblem((y, p, x) -> dy_dx(x), y1, (x1, x2)), Tsit5(); reltol=1e-10)
    xs, ys = sol.t, sol.u
    return Spline1D(xs, ys; k=3, bc="error"), x1, x2
end

# conformal time
function η(co::ΛCDM, x::Real)
    if isnothing(co.η_spline)
        dη_dx(x) = 1 / (a(x) * H(co, x))
        x1, x2 = -20.0, +20.0
        aeq = co.Ωr0 / co.Ωm0
        if co.Ωm0 > 0
            η1 = 2 / (co.H0*√(co.Ωm0)) * (√(a(x1)+aeq) - √(aeq))
        else
            η1 = 1 / aH(co, x1)
        end
        co.η_spline, x1, x2 = _spline_dy_dx(co, dη_dx, x1, x2, η1)
    end
    (x1, x2) = extrema(get_knots(co.η_spline))
    return x1 <= x <= x2 ? co.η_spline(x) : NaN
end

# cosmic time
function t(co::ΛCDM, x::Real)
    if isnothing(co.t_spline)
        dt_dx(x) = 1 / H(co, x)
        x1, x2 = -20.0, +20.0 # integration and spline range
        aeq = co.Ωr0 / co.Ωm0
        if co.Ωm0 > 0
            t1 = 2 / (3*co.H0*√(co.Ωm0)) * (√(a(x1)+aeq) * (a(x1)-2*aeq) + 2*aeq^(3/2))
        else
            t1 = 1 / (2*H(co, x1))
        end
        co.t_spline, x1, x2 = _spline_dy_dx(co, dt_dx, x1, x2, t1)
    end
    (x1, x2) = extrema(get_knots(co.t_spline))
    return x1 <= x <= x2 ? co.t_spline(x) : NaN
end

# density parameters (relative to critical density *at the time*)
# computed using Ωs = ρs/ρcrit = ρs/ρcrit0 * ρcrit0/ρcrit = Ωs0 * H0^2/H^2
Ωγ(co::ΛCDM, x::Real) = co.Ωγ0 / (a(x)^4 * H(co, x)^2 / H(co, 0)^2)
Ων(co::ΛCDM, x::Real) = co.Ων0 / (a(x)^4 * H(co, x)^2 / H(co, 0)^2)
Ωb(co::ΛCDM, x::Real) = co.Ωb0 / (a(x)^3 * H(co, x)^2 / H(co, 0)^2)
Ωc(co::ΛCDM, x::Real) = co.Ωc0 / (a(x)^3 * H(co, x)^2 / H(co, 0)^2)
Ωk(co::ΛCDM, x::Real) = co.Ωk0 / (a(x)^2 * H(co, x)^2 / H(co, 0)^2)
ΩΛ(co::ΛCDM, x::Real) = co.ΩΛ  / (H(co, x)^2 / H(co, 0)^2)
Ωr(co::ΛCDM, x::Real) = Ωγ(co, x) + Ων(co, x)
Ωm(co::ΛCDM, x::Real) = Ωb(co, x) + Ωc(co, x)
Ω( co::ΛCDM, x::Real) = Ωr(co, x) + Ωm(co, x) + Ωk(co, x) + ΩΛ(co, x)

# equalities
r_m_equality(co::ΛCDM) = log(co.Ωr0 / co.Ωm0)
m_Λ_equality(co::ΛCDM) = log(co.Ωm0 / co.ΩΛ) / 3
acceleration_onset(co::ΛCDM) = find_zero(x -> daH(co, x), (-20, +20))

# luminosity distance
function dL(co::ΛCDM, x::Real)
    χ = c * (η(co, 0) - η(co, x))
    # 1. r(Ωk0 = 0) = χ
    # 2. r(Ωk0 < 0) = χ *  sin(√(-Ωk0)*H0*χ/c) / (√(-Ωk0)*H0*χ/c)
    # 3. r(Ωk0 > 0) = χ * sinh(√(+Ωk0)*H0*χ/c) / (√(+Ωk0)*H0*χ/c)
    # can all be written as the real-valued function
    #    r(Ωk0)     = χ * sinc(√(-Ωk0)*H0*χ/c/π)
    # using complex numbers,
    # because sinc(x) = sin(π*x) / (π*x) -> 1 as x -> 0,
    # and sinh(x) = -i * sin(i*x)
    r = χ * real(sinc(√(complex(-co.Ωk0)) * co.H0 * χ / c / π)) # sinc(x) = sin(π*x) / (π*x), so divide argument by π!
    return r / a(x)
end

function dA(co::ΛCDM, x::Real)
    return dL(co, x) * a(x)^2
end

# EXAMPLES:
# Cosmology.ΛCDM(Ωr0=0, Ωb0=0, Ωc0=0.2, Ωk0=-0.9) has Ωpoly(x ≈ -1) < 0
function is_fucked(co::ΛCDM; x1=-20.0, x2=+20.0)
    if Ωpoly(co, x1) < 0 || Ωpoly(co, x2) < 0
        return true
    end

    a = -4 * co.Ωr0
    b = -3 * co.Ωm0
    c = -2 * co.Ωk0
    d = b^2 - 4*a*c
    if d >= 0
        ainv1 = (-b + √(d)) / (2*a)
        ainv2 = (-b - √(d)) / (2*a)
        a1, a2 = 1/ainv1, 1/ainv2
        if (a1 >= 0 && Ωpoly(co, x(a1)) < 0) || (a2 >= 0 && Ωpoly(co, x(a2)) < 0)
            return true
        end
    end

    return false
end

end
