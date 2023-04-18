module Cosmology

include("Constants.jl")

using .Constants
using Roots # root finding
#using OrdinaryDiffEq # ODE integration (instead of DifferentialEquations to reduce compile time)
using DifferentialEquations

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

export time_tight_coupling
export δc, δb, vc, vb, Φ, Θl

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
          x_switch_Peebles::Float64

    # splines (lazily initialized)
    η_spline::Union{Nothing, SciMLBase.ODESolution} # conformal time
    t_spline::Union{Nothing, SciMLBase.ODESolution} # cosmic    time
    Xe_Peebles_spline::Union{Nothing, SciMLBase.ODESolution} # free electron fraction (TODO: separate struct?)
    τ_spline::Union{Nothing, SciMLBase.ODESolution} # optical depth (TODO: separate struct?)
    perturbations_tight_spline::Union{Nothing, SciMLBase.ODESolution} # perturbations during tight coupling
    perturbations_untight_spline::Union{Nothing, SciMLBase.ODESolution} # perturbations after tight coupling

    function ΛCDM(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046, Yp=0.24, z_reion_H=8.0, Δz_reion_H=0.5, z_reion_He=3.5, Δz_reion_He=0.5)
        H0  = h * 100*km/Mpc # 1/s
        Ωm0 = Ωb0 + Ωc0
        Ωγ0 = π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ0 = 1.0 - (Ωr0 + Ωm0 + Ωk0)
        reionization = count(isnan(num) for num in (z_reion_H, z_reion_He, Δz_reion_H, Δz_reion_He)) == 0 # turn off by setting either to NaN
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ0, Tγ0, Neff, Yp, reionization, z_reion_H, Δz_reion_H, z_reion_He, Δz_reion_He, NaN, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end

# in vectorized calls, like H.(co, [1.0, 2.0]),
# broadcast the same cosmology to all scalar calls
Base.broadcastable(co::ΛCDM) = Ref(co)

include("Background.jl")
include("Recombination.jl")
include("Perturbations.jl")
include("Utils.jl")

end
