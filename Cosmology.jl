module Cosmology

include("Constants.jl")

using .Constants # physical constants
using Roots # root finding
using DifferentialEquations # ODE integration
using ODEInterfaceDiffEq # fast & stiff "radau" ODE integrator (for perturbations)
using Dierckx # splines
using Base.Threads # parallelization
using SpecialFunctions: sphericalbesselj as jl # spherical Bessel function # TODO: correct function?

export a, z
export ΛCDM
export H, aH, daH, d2aH
export t, η, Ωγ, Ων, Ωb, Ωc, Ωk, ΩΛ, Ωr, Ωm, Ω
export equality_rm, equality_mΛ, acceleration_onset
export dL, dA

export Tγ, Xe_Saha_H, Xe_Saha_H_He, Xe_Peebles, Xe, τ, dτ, d2τ, g, dg, d2g, time_switch_Peebles
export time_decoupling, time_recombination, sound_horizon
export time_reionization_H, time_reionization_He
export multirange
export format_time_variations

export time_tight_coupling, time_horizon_entry
export δc, δb, vc, vb, Φ, Ψ, Θl, Nl, ΘPl, S

export P_primordial, P

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

    # recombination # TODO: separate struct?
    const Yp::Float64 # helium mass fraction
    const reionization::Bool   # reionization on?
    const  z_reion_H::Float64  # Hydrogen reionization redshift time
    const Δz_reion_H::Float64  # Hydrogen reionization redshift duration
    const  z_reion_He::Float64 # Helium   reionization redshift time
    const Δz_reion_He::Float64 # Helium   reionization redshift duration

    # primordial power spectrum
    const As::Float64 # power spectrum amplitude (for k = k_pivot)
    const ns::Float64 # power spectrum amplitude (for k = k_pivot)
    const k_pivot::Float64

    # splines (lazily initialized)
    η_spline::Union{Nothing, Spline1D} # conformal time
    t_spline::Union{Nothing, Spline1D} # cosmic    time
    Xe_spline::Union{Nothing, Spline1D} # free electron fraction (TODO: separate struct?)
    τ_spline::Union{Nothing, Spline1D} # optical depth (TODO: separate struct?)
    g_spline::Union{Nothing, Spline1D} # optical depth (TODO: separate struct?)
    sound_horizon_spline::Union{Nothing, Spline1D} # sound horizon

    perturbation_splines1D::Vector{Tuple{Float64, Vector{Spline1D}}} # cached (k, [y1(x), y2(x), ...]) pairs
    perturbation_splines2D::Vector{Union{Nothing, Spline2D}} # [y1(x, k), y2(x, k), ...]

    function ΛCDM(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046, Yp=0.24, z_reion_H=8.0, Δz_reion_H=0.5, z_reion_He=3.5, Δz_reion_He=0.5, As=2e-9, ns=0.96, k_pivot=0.05/Mpc)
        H0  = h * 100*km/Mpc # 1/s
        Ωm0 = Ωb0 + Ωc0
        Ωγ0 = π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ0 = 1.0 - (Ωr0 + Ωm0 + Ωk0)
        reionization = count(isnan(num) for num in (z_reion_H, z_reion_He, Δz_reion_H, Δz_reion_He)) == 0 # turn off by setting either to NaN
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ0, Tγ0, Neff, Yp, reionization, z_reion_H, Δz_reion_H, z_reion_He, Δz_reion_He, As, ns, k_pivot, nothing, nothing, nothing, nothing, nothing, nothing, [], [])
    end
end

# in vectorized calls, like H.(co, [1.0, 2.0]),
# broadcast the same cosmology to all scalar calls
Base.broadcastable(co::ΛCDM) = Ref(co)

include("Background.jl")
include("Recombination.jl")
include("Perturbations.jl")
include("Observables.jl")
include("Utils.jl")

end
