include("Constants.jl")

using OrdinaryDiffEq # ODE integration (instead of DifferentialEquations to reduce compile time)
using Dierckx # spline interpolation
using Roots # root finding
using .Constants

mutable struct ΛCDM
    const h0::Float64  # dimensionless Hubble parameter h0 = H0 / (100 (km/s/Mpc))
    const H0::Float64  # dimensionful  Hubble parameter
    const Ωb0::Float64 # baryons
    const Ωc0::Float64 # cold dark matter
    const Ωm0::Float64 # matter = baryons + cold dark matter
    const Ωk0::Float64 # curvature
    const Ωγ0::Float64 # photons
    const Ων0::Float64 # neutrinos (massless)
    const Ωr0::Float64 # radiation = photons + neutrinos
    const ΩΛ0::Float64 # dark energy (as cosmological constant)
    const Tγ0::Float64 # photon (CMB) temperature
    const Neff::Float64 # effective neutrino number

    η_spline::Union{Nothing, Dierckx.Spline1D} # conformal time spline (lazily initialized)
    t_spline::Union{Nothing, Dierckx.Spline1D} # cosmic    time spline (lazily initialized)

    function ΛCDM(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046)
        H0  = h * 100*km/Mpc # 1/s
        Ωm0 = Ωb0 + Ωc0
        Ωγ0 = π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ0 = 1.0 - (Ωr0 + Ωm0 + Ωk0)
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ0, Tγ0, Neff, nothing, nothing)
    end
end

# in vectorized calls, like H.(co, [1.0, 2.0]),
# broadcast the same cosmology to all scalar calls
Base.broadcastable(co::ΛCDM) = Ref(co)

# convert between (a, x, z) across a = e^x = 1/(z+1)
x(a::Real) = log(a)     # natural logarithm of scale factor (our internal time variable)
a(x::Real) = exp(x)     # scale factor
z(x::Real) = 1/a(x) - 1 # redshift

# the d'th derivative of the thing in the square root in the Friedmann equation, E = H^2 / H0^2
E(co::ΛCDM, x::Real; d::Integer=0) = (-4)^d * co.Ωr0/a(x)^4 +
                                     (-3)^d * co.Ωm0/a(x)^3 +
                                     (-2)^d * co.Ωk0/a(x)^2 +
                                        0^d * co.ΩΛ0          # 0^0 = 1, 0^1 = 0, 0^2 = 0, ...

# Friedmann equation
   H(co::ΛCDM, x::Real) = co.H0 * √(E(co, x)) # cosmic    Hubble parameter
  aH(co::ΛCDM, x::Real) = a(x) * H(co, x)     # conformal Hubble parameter
 daH(co::ΛCDM, x::Real) = aH(co, x) * (1 + 1/2 * E(co, x; d=1) / E(co, x)) # 1st derivative of aH
d2aH(co::ΛCDM, x::Real) = aH(co, x) * (1 + E(co, x; d=1) / E(co, x) + 1/2 * E(co, x; d=2) / E(co, x) - 1/4 * (E(co, x; d=1) / E(co, x))^2) # 2nd derivative of aH

# internal function that computes y(x) on a cubic spline,
# given dy/dx, from x=x1 with y(x1)=y1, to x=x2
function _spline_dy_dx(co::ΛCDM, dy_dx::Function, x1::Float64, x2::Float64, y1::Float64)
    sol = solve(ODEProblem((y, p, x) -> dy_dx(x), y1, (x1, x2)), Tsit5(); reltol=1e-10)
    xs, ys = sol.t, sol.u
    return Spline1D(xs, ys; k=3, bc="error") # throw error if evaluating spline outside bounds
end

# conformal time
function η(co::ΛCDM, x::Real)
    if isnothing(co.η_spline)
        # lazy initialize spline
        dη_dx(x) = 1 / aH(co, x)
        x1, x2 = -20.0, +20.0
        aeq = co.Ωr0 / co.Ωm0
        if co.Ωm0 > 0
            η1 = 2 / (co.H0*√(co.Ωm0)) * (√(a(x1)+aeq) - √(aeq)) # anal expr with Ωk=ΩΛ=0
        else
            η1 = 1 / aH(co, x1) # anal expr with Ωm=Ωk=ΩΛ=0
        end
        co.η_spline = _spline_dy_dx(co, dη_dx, x1, x2, η1)
    end
    return co.η_spline(x)
end

# cosmic time
function t(co::ΛCDM, x::Real)
    if isnothing(co.t_spline)
        # lazy initialize spline
        dt_dx(x) = 1 / H(co, x)
        x1, x2 = -20.0, +20.0
        aeq = co.Ωr0 / co.Ωm0
        if co.Ωm0 > 0
            t1 = 2 / (3*co.H0*√(co.Ωm0)) * (√(a(x1)+aeq) * (a(x1)-2*aeq) + 2*aeq^(3/2)) # anal expr with Ωk=ΩΛ=0
        else
            t1 = 1 / (2*H(co, x1)) # anal expr with Ωm=Ωk=ΩΛ=0
        end
        co.t_spline = _spline_dy_dx(co, dt_dx, x1, x2, t1)
    end
    return co.t_spline(x)
end

# density parameters (relative to critical density *at the time*)
# computed using Ωs = ρs/ρcrit = ρs/ρcrit0 * ρcrit0/ρcrit = Ωs0 * H0^2/H^2 = Ωs0 / E^2
Ωγ(co::ΛCDM, x::Real) = co.Ωγ0 / a(x)^4 / E(co, x)
Ων(co::ΛCDM, x::Real) = co.Ων0 / a(x)^4 / E(co, x)
Ωb(co::ΛCDM, x::Real) = co.Ωb0 / a(x)^3 / E(co, x)
Ωc(co::ΛCDM, x::Real) = co.Ωc0 / a(x)^3 / E(co, x)
Ωk(co::ΛCDM, x::Real) = co.Ωk0 / a(x)^2 / E(co, x)
ΩΛ(co::ΛCDM, x::Real) = co.ΩΛ0          / E(co, x)
Ωr(co::ΛCDM, x::Real) = Ωγ(co, x) + Ων(co, x)
Ωm(co::ΛCDM, x::Real) = Ωb(co, x) + Ωc(co, x)
Ω( co::ΛCDM, x::Real) = Ωr(co, x) + Ωm(co, x) + Ωk(co, x) + ΩΛ(co, x)

# time of equality between different species (as x = log(a))
equality_rm(co::ΛCDM) = log(co.Ωr0 / co.Ωm0)
equality_mΛ(co::ΛCDM) = log(co.Ωm0 / co.ΩΛ0) / 3

# time of acceleration onset (as x = log(a))
acceleration_onset(co::ΛCDM) = find_zero(x -> daH(co, x), (-20, +20))

# conformal distance
χ(co::ΛCDM, x::Real) = c * (η(co, 0) - η(co, x))

# radial coordinate (of light emitted at x)
#   note: ALL three expressions
#   1. r(Ωk0 = 0) = χ
#   2. r(Ωk0 < 0) = χ *  sin(√(-Ωk0)*H0*χ/c) / (√(-Ωk0)*H0*χ/c)
#   3. r(Ωk0 > 0) = χ * sinh(√(+Ωk0)*H0*χ/c) / (√(+Ωk0)*H0*χ/c)
#   can be written as the real-valued function
#      r(Ωk0)     = χ * sinc(√(-Ωk0)*H0*χ/c/π)
#   using complex numbers,
#   because sinc(x) = sin(π*x) / (π*x) -> 1 as x -> 0,
#   and sinh(x) = -i * sin(i*x)
r(co::ΛCDM, x::Real) = χ(co, x) * real(sinc(√(complex(-co.Ωk0)) * co.H0 * χ(co, x) / c / π)) # in Julia, sinc(x) = sin(π*x) / (π*x), so divide argument by π!

# angular diameter distance and luminosity distance (of light emitted at x)
dA(co::ΛCDM, x::Real) = r(co, x) * a(x)
dL(co::ΛCDM, x::Real) = r(co, x) / a(x)

# checks whether the Hubble parameter becomes complex on the integration interval (x1, x2)
# EXAMPLE: Cosmology.ΛCDM(Ωr0=0, Ωb0=0, Ωc0=0.2, Ωk0=-0.9) has E(x ≈ -1) < 0
function is_fucked(co::ΛCDM; x1=-20.0, x2=+20.0)
    # negative at the endpoints?
    if E(co, x1) < 0 || E(co, x2) < 0
        return true
    end

    # negative between the endpoints?
    # check the values at the stationary points of E(x), defined by roots of 2nd degree polynomial
    # dE/dx = -a^(-2) * [ 4*Ωr0*a^(-2) + 3*Ωm0*a^(-1) + 2*Ωk0 ] = 0
    a = 4 * co.Ωr0
    b = 3 * co.Ωm0
    c = 2 * co.Ωk0
    d = b^2 - 4*a*c
    if d >= 0
        ainv1 = (-b + √(d)) / (2*a)
        ainv2 = (-b - √(d)) / (2*a)
        a1, a2 = 1/ainv1, 1/ainv2
        if (a1 >= 0 && E(co, x(a1)) < 0) || (a2 >= 0 && E(co, x(a2)) < 0)
            return true
        end
    end

    return false # always positive
end
