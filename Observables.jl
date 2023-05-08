#=
struct PowerSpectrum
    const As::Float64 # power spectrum amplitude (for k = k_pivot)
    const ns::Float64 # power spectrum spectral index (for k = k_pivot)
    const k_pivot::Float64

    function PowerSpectrum(; As=2e-9, ns=0.96, k_pivot=0.05/Mpc)
        new(As, ns, k_pivot)
    end
end

P_primordial(co::ΛCDM, k::Real) = 2*π^2/k^3 * co.As * (k/co.k_pivot)^(co.ns-1)
Δ(co::ΛCDM, x::Real, k::Real) = c^2*k^2*Φ(co,x,k) / (3/2*co.Ωm0/a(x)*co.H0^2)
P(co::ΛCDM, x::Real, k::Real) = abs(Δ(co,x,k))^2 * P_primordial(co, k)

function Cl(co::ΛCDM, l::Int)
    println("l = $l")
    kmin =    1 / (c*η(co,0))
    kmax = 3000 / (c*η(co,0))
    Δk = 2*π / (c*η(co,0)*10)
    k = range(kmin, kmax, step=Δk)
    x = range(-20, 0, length=1000)
    Θl0(l, k) = trapz(x, S(co,x,k) .* jl.(l, c*k*(η(co,0.0) .- η.(co,x)))) # TODO: or trapz?
    return 2/π * trapz(k, @. k^2 * P_primordial(co,k) * Θl0(l,k)^2)
end
=#
