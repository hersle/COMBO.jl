P_primordial(co::ΛCDM, k::Real) = 2*π^2/k^3 * co.As * (k/co.k_pivot)^(co.ns-1)
Δ(co::ΛCDM, x::Real, k::Real) = c^2*k^2*Φ(co,x,k) / (3/2*co.Ωm0/a(x)*co.H0^2)
P(co::ΛCDM, x::Real, k::Real) = abs(Δ(co,x,k))^2 * P_primordial(co, k)
