include("Constants.jl")

ρcrit(co::ΛCDM, x::Real) = 3 * H(co, x)^2 / (8 * π * G)
ρb(co::ΛCDM, x::Real) = Ωb(co, x) * ρcrit(co, x)
nH(co::ΛCDM, x::Real; Yp::Real=0) = (1-Yp) * ρb(co, x) / mH # TODO: need to change when adding He?
Tγ(co::ΛCDM, x::Real) = co.Tγ0 / a(x)

function Xe_Saha(co::ΛCDM, x::Real)
    # TODO: add He
    # TODO: fix overflow issues!
    T = Tγ(co, x)
    a = 1
    nb = nH(co,x) # TODO: generalize to helium
    b = (2 * π * me * kB * T / h^2)^(3/2) * exp(-EHion/(kB*T)) / nb
    c = -b
    return quadroots(a, b, c)[2]
end

function Xe_Peebles(co::ΛCDM, x::Real, x1::Real, Xe1::Real)
    if isnothing(co.Xe_Peebles_spline)
        function dXe_dx(x, Xe)
            Yp = 0.0 # TODO: add He
            Tb = Tγ(co,x) # K # TODO: assumptions?
            n_1s = (1-Xe) * nH(co,x;Yp) # 1/m^3
            Λ_2s_1s = 8.227 # 1/s
            Λ_α = H(co,x) * (3*EHion/(ħ*c))^3 / ((8*π)^2 * n_1s) # 1/s
            ϕ2 = 0.448 * log(EHion/(kB*Tb)) # dimensionless
            α2 = 64*π / √(27*π) * (α/me)^2 * √(EHion/(kB*Tb)) * ϕ2 * ħ^2/c # m^3/s
            λe = √(h^2 / (2*π*me*kB*Tb)) # thermal de Broglie wavelength
            β  = α2 / λe^3 * exp(-EHion/(kB*Tb))
            β2 = α2 / λe^3 * exp(-EHion/(4*kB*Tb)) # 1/s (compute this instead of β2 = β * exp(3*EHion/(4*kB*Tb)) to avoid exp overflow)
            C_r = (Λ_2s_1s + Λ_α) / (Λ_2s_1s + Λ_α + β2)
            return C_r / H(co,x) * (β*(1-Xe) - nH(co,x;Yp)*α2*Xe^2)
        end
        co.Xe_Peebles_spline = _spline_integral(dXe_dx, x1, +20.0, Xe1)
    end

    return co.Xe_Peebles_spline(x) # TODO: spline the logarithm instead?
end

function Xe(co::ΛCDM, x::Real)
    Xe1 = 0.99
    x1  = find_zero(x -> Xe_Saha(co, x) - Xe1, (-20.0, +20.0))

    if x < x1
        return Xe_Saha(co, x) # regime where Saha equation is valid
    else
        return Xe_Peebles(co, x, x1, Xe1) # regime where Peebles equation is valid
    end
end

ne(co::ΛCDM, x::Real) = nH(co,x) * Xe(co,x)

 dτ(co::ΛCDM, x::Real) = -ne(co,x) * σT * c / H(co,x)
d2τ(co::ΛCDM, x::Real; Δx::Real=1e-5) = (dτ(x+Δx/2) - dτ(x-Δx/2)) / Δx

function τ(co::ΛCDM, x::Real)
    if isnothing(co.τ_spline)
        co.τ_spline = _spline_integral((x, τ) -> dτ(co, x), 0.0, -20.0, 0.0)
    end
    return co.τ_spline(x)
end
