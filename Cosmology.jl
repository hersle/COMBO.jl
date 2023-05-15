module Cosmology

# TODO: use DifferentialEquations' splines
# TODO: use autodiff

include("Constants.jl")
using .Constants # physical constants

using Roots # root finding
using DifferentialEquations # ODE integration
using ForwardDiff # compile-time analytical derivatives
using Base.Threads # parallelization
using Bessels: sphericalbesselj as jl # spherical Bessel function # TODO: correct function?
using Trapz, QuadGK # quadrature
using Interpolations # splines
using OffsetArrays # arbitrary-indexed vectors

export a, z
export Parameters, Background, Recombination, PerturbationMode
export H, aH, aH′, aH′′
export t, η, Ωγ, Ων, Ωb, Ωc, Ωk, ΩΛ, Ωr, Ωm, Ω
export equality_rm, equality_mΛ, acceleration_onset
export dL, dA

export Tγ, Xe_Saha_H, Xe_Saha_H_He, Xe_Peebles, Xe, τ, τ′, τ′′, g, g′, g′′, time_switch_Peebles
export time_decoupling, time_recombination, s
export time_reionization_H, time_reionization_He
export multirange
export format_time_variations

export time_tight_coupling, time_horizon_entry
export δc, δb, vc, vb, Φ, Ψ, Θl, Nl, ΘPl, S

export P_primordial, P, Cl

struct Parameters
    # background parameters
    h0::Float64  # dimensionless Hubble parameter h0 = H0 / (100 (km/s/Mpc))
    H0::Float64  # dimensionful  Hubble parameter
    Ωb0::Float64 # baryons
    Ωc0::Float64 # cold dark matter
    Ωm0::Float64 # matter = baryons + cold dark matter
    Ωk0::Float64 # curvature
    Ωγ0::Float64 # photons
    Ων0::Float64 # neutrinos (massless)
    Ωr0::Float64 # radiation = photons + neutrinos
    ΩΛ0::Float64 # dark energy (as cosmological constant)
    Tγ0::Float64 # photon (CMB) temperature
    Neff::Float64 # effective neutrino number

    # recombination parameters
    Yp::Float64 # helium mass fraction
    reionization::Bool   # reionization on?
     z_reion_H::Float64  # Hydrogen reionization redshift time
    Δz_reion_H::Float64  # Hydrogen reionization redshift duration
     z_reion_He::Float64 # Helium   reionization redshift time
    Δz_reion_He::Float64 # Helium   reionization redshift duration

    # power spectrum parameters
    As::Float64 # power spectrum amplitude (for k = k_pivot)
    ns::Float64 # power spectrum spectral index (for k = k_pivot)
    k_pivot::Float64

    function Parameters(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046, Yp=0.24, z_reion_H=8.0, Δz_reion_H=0.5, z_reion_He=3.5, Δz_reion_He=0.5, As=2e-9, ns=0.96, k_pivot=0.05/Mpc)
        H0  = h * 100*km/Mpc # 1/s
        Ωm0 = Ωb0 + Ωc0
        Ωγ0 = π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ0 = 1.0 - (Ωr0 + Ωm0 + Ωk0)
        reionization = all(!isnan(num) for num in (z_reion_H, z_reion_He, Δz_reion_H, Δz_reion_He)) # turn off by setting either to NaN
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ0, Tγ0, Neff, Yp, reionization, z_reion_H, Δz_reion_H, z_reion_He, Δz_reion_He, As, ns, k_pivot)
    end
end

Base.broadcastable(par::Parameters) = Ref(par)

# convert between (a, x, z) across a = e^x = 1/(z+1)
x(a) = log(a)     # natural logarithm of scale factor (our internal time variable)
a(x) = exp(x)     # scale factor
z(x) = 1/a(x) - 1 # redshift

# the d'th derivative of the thing in the square root in the Friedmann equation, E = H^2 / H0^2
# TODO: does this also work with ForwardDiff?
E(par::Parameters, x; d=0) = (-4)^d * par.Ωr0/a(x)^4 +
                             (-3)^d * par.Ωm0/a(x)^3 +
                             (-2)^d * par.Ωk0/a(x)^2 +
                             0^d * par.ΩΛ0          # 0^0 = 1, 0^1 = 0, 0^2 = 0, ...

# Friedmann equation
H(par::Parameters, x) = par.H0 * √(E(par, x)) # cosmic    Hubble parameter
aH(par::Parameters, x) = a(x) * H(par, x)     # conformal Hubble parameter
aH′(par::Parameters, x) = ForwardDiff.derivative(x ->  aH(par, x), x)
aH′′(par::Parameters, x) = ForwardDiff.derivative(x -> aH′(par, x), x)

# density parameters (relative to critical density *at the time*)
# computed using Ωs = ρs/ρcrit = ρs/ρcrit0 * ρcrit0/ρcrit = Ωs0 * H0^2/H^2 = Ωs0 / E
Ωγ(par::Parameters, x) = par.Ωγ0 / a(x)^4 / E(par, x)
Ων(par::Parameters, x) = par.Ων0 / a(x)^4 / E(par, x)
Ωb(par::Parameters, x) = par.Ωb0 / a(x)^3 / E(par, x)
Ωc(par::Parameters, x) = par.Ωc0 / a(x)^3 / E(par, x)
Ωk(par::Parameters, x) = par.Ωk0 / a(x)^2 / E(par, x)
ΩΛ(par::Parameters, x) = par.ΩΛ0          / E(par, x)
Ωr(par::Parameters, x) = Ωγ(par, x) + Ων(par, x)
Ωm(par::Parameters, x) = Ωb(par, x) + Ωc(par, x)
Ω( par::Parameters, x) = Ωr(par, x) + Ωm(par, x) + Ωk(par, x) + ΩΛ(par, x)

aeq(par::Parameters) = par.Ωr0 / par.Ωm0

# time of equality between different species (as x = log(a))
# TODO: rename time_...
equality_rm(par::Parameters) = log(par.Ωr0 / par.Ωm0)
equality_mΛ(par::Parameters) = log(par.Ωm0 / par.ΩΛ0) / 3

# time of acceleration onset (as x = log(a))
acceleration_onset(par::Parameters) = find_zero(x -> aH′(par, x), (-20, +20))

# checks whether the Hubble parameter becomes complex on the integration interval (x1, x2)
# EXAMPLE: Cosmology.ΛCDM(Ωr0=0, Ωb0=0, Ωc0=0.2, Ωk0=-0.9) has E(x ≈ -1) < 0
function is_fucked(par::Parameters; x1=-20.0, x2=+20.0)
    # negative at the endpoints?
    if E(par, x1) < 0 || E(par, x2) < 0
        return true
    end

    # negative between the endpoints?
    # check the values at the stationary points of E(x), defined by roots of 2nd degree polynomial
    # dE/dx = -a^(-2) * [ 4*Ωr0*a^(-2) + 3*Ωm0*a^(-1) + 2*Ωk0 ] = 0
    ainv1, ainv2 = quadroots(4*par.Ωr0, 3*par.Ωm0, 2*par.Ωk0)
    if !isnan(ainv1)
        a1, a2 = 1/ainv1, 1/ainv2
        if (a1 >= 0 && E(par, x(a1)) < 0) || (a2 >= 0 && E(par, x(a2)) < 0)
            return true
        end
    end

    return false # always positive
end

include("Background.jl")
include("Recombination.jl")
include("Perturbations.jl")
include("Observables.jl")
include("Utils.jl")

end
