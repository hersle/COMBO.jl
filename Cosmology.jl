module Cosmology

include("Constants.jl")
using .Constants # physical constants

using Roots # root finding
using OrdinaryDiffEq # ODE integration (part of DifferentialEquations)
using ForwardDiff # compile-time analytical derivatives

export Parameters
export a, z

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
    fHe::Float64 # helium-to-hydrogen ratio
    reionization::Bool   # reionization on?
     z_reion_H::Float64  # Hydrogen reionization redshift time
    Δz_reion_H::Float64  # Hydrogen reionization redshift duration
     z_reion_He::Float64 # Helium   reionization redshift time
    Δz_reion_He::Float64 # Helium   reionization redshift duration

    # power spectrum parameters
    As::Float64 # power spectrum amplitude (for k = k_pivot)
    ns::Float64 # power spectrum spectral index (for k = k_pivot)
    k_pivot::Float64

    function Parameters(; h=0.67, Ωb0=0.05, Ωc0=0.267, Ωk0=0, Tγ0=2.7255, Neff=3.046, Yp=0.245, z_reion_H=8.0, Δz_reion_H=0.5, z_reion_He=3.5, Δz_reion_He=0.5, As=2.1e-9, ns=0.965, k_pivot=0.05/Mpc)
        H0  = h * 100*km/Mpc # 1/s
        Ωm0 = Ωb0 + Ωc0
        Ωγ0 = π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8*π*G / (3*H0^2)
        Ων0 = Neff * 7/8 * (4/11)^(4/3) * Ωγ0
        Ωr0 = Ωγ0 + Ων0
        ΩΛ0 = 1.0 - (Ωr0 + Ωm0 + Ωk0)
        fHe = Yp / (4*(1-Yp))
        reionization = all(!isnan(num) for num in (z_reion_H, z_reion_He, Δz_reion_H, Δz_reion_He)) # turn off by setting either to NaN
        new(h, H0, Ωb0, Ωc0, Ωm0, Ωk0, Ωγ0, Ων0, Ωr0, ΩΛ0, Tγ0, Neff, Yp, fHe, reionization, z_reion_H, Δz_reion_H, z_reion_He, Δz_reion_He, As, ns, k_pivot)
    end
end

# convert between (a, x, z) across a = e^x = 1/(z+1)
x(a) = log(a)     # natural logarithm of scale factor (our internal time variable)
a(x) = exp(x)     # scale factor
z(x) = 1/a(x) - 1 # redshift

include("Background.jl")
include("Recombination.jl")
include("Perturbations.jl")
include("PowerSpectra.jl")
include("Utils.jl")

# keep structs fixed in vectorized calls like aH.(par, xs)
Base.broadcastable(s::Union{Parameters,Background,Recombination,PerturbationMode}) = Ref(s)

end
