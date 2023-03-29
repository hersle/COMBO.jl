include("Constants.jl")

ρcrit(co::ΛCDM, x::Real) = 3 * H(co, x)^2 / (8 * π * G)
ρb(co::ΛCDM, x::Real) = Ωb(co, x) * ρcrit(co, x)
nb(co::ΛCDM, x::Real) = ρb(co, x) / mH
nH(co::ΛCDM, x::Real) = (1-co.Yp) * nb(co, x)
nHe(co::ΛCDM, x::Real) = co.Yp/4 * nb(co, x)
Tγ(co::ΛCDM, x::Real) = co.Tγ0 / a(x)
fHe(Yp::Real) = Yp / (4*(1-Yp)) # TODO: what is this? replace with my hlines?

function Xe_Saha_H(co::ΛCDM, x::Real)
    @assert co.Yp == 0

    T = Tγ(co, x)
    a = 1
    nb = nH(co,x)
    b = (2 * π * me * kB * T / h^2)^(3/2) * exp(-EH1ion/(kB*T)) / nb
    c = -b
    #return Float64(quadroots(BigFloat(a), BigFloat(b), BigFloat(c))[2])
    # when b >> 1, the quadratic equation solution is
    #   (-b + √(b^2+4*b)) / 2
    # = b/2 * (-1 + √(1+4/b))
    # ≈ b/2 * (-1 + 1 + 2/b)  (b >> 1)
    # = 1
    return b < 1e10 ? quadroots(a, b, c)[2] : 1.0
end

# TODO: slow/fails for x ≲ -7
function Xe_Saha_H_He(co::ΛCDM, x::Real; tol::Float64=1e-15, maxiters::Int=10000)
    # To find Xe = XH+ + Yp/(4*(1-Yp)) * (XHe+ + 2*XHe++),
    # begin with an initial guess for Xe and iteratively solve the system of Saha equations
    # (1) ne * XHe+ / (1 - XHe+ - XHe++) = 2 / λe^3 * exp(-EHe1ion / (kB*Tb))
    # (2) ne * XHe++ / XHe+ = 4 / λe^3 * exp(-EHe2ion/(kB*Tb))
    # (3) ne * XH+ / (1-XH+) = 1 / λe^3 * exp(-EH1ion/(kB*Tb))
    # for {XH+, XHe+, XHe++} = {XH1, XHe1, XHe2}.

    Tb = Tγ(co, x)
    λe = h / √(2*π*me*kB*Tb)
    RHS1 = 2 / λe^3 * exp(-EHe1ion/(kB*Tb))
    RHS2 = 4 / λe^3 * exp(-EHe2ion/(kB*Tb))
    RHS3 = 1 / λe^3 * exp(-EH1ion /(kB*Tb))

    # 1. Guess Xe
    Xe = 1.0
    converged = false
    iterations = 0
    while !converged && iterations < maxiters
        # 2. Compute corresponding XH+, XHe+, XHe++
        ne = Xe * nH(co, x)
        XH1  = RHS3 / (ne + RHS3)
        XHe1 = 1 / (1 + ne/RHS1 + RHS2/ne)
        XHe2 = RHS2 / ne * XHe1

        # 3. Compute corresponding Xe ...
        Xe_new = XH1 + fHe(co.Yp) * (XHe1 + 2*XHe2)
        converged = abs(Xe_new - Xe) < tol # ... until Xe becomes self consistent
        Xe = Xe_new
        iterations += 1
    end

    # println("$iterations iterations")
    return iterations < maxiters ? Xe : NaN
end

function Xe_Peebles(co::ΛCDM, x::Real, x1::Real, Xe1::Real)
    if isnothing(co.Xe_Peebles_spline)
        function dXe_dx(x, Xe)
            Tb = Tγ(co,x) # K # TODO: assumptions? separate baryon evolution?
            n_1s = (1-Xe) * nH(co,x) # 1/m^3
            Λ_2s_1s = 8.227 # 1/s
            Λ_α = H(co,x) * (3*EH1ion/(ħ*c))^3 / ((8*π)^2 * n_1s) # 1/s
            ϕ2 = 0.448 * log(EH1ion/(kB*Tb)) # dimensionless
            α2 = 64*π / √(27*π) * (α/me)^2 * √(EH1ion/(kB*Tb)) * ϕ2 * ħ^2/c # m^3/s
            λe = √(h^2 / (2*π*me*kB*Tb)) # thermal de Broglie wavelength
            β  = α2 / λe^3 * exp(-EH1ion/(kB*Tb))
            β2 = α2 / λe^3 * exp(-EH1ion/(4*kB*Tb)) # 1/s (compute this instead of β2 = β * exp(3*EH1ion/(4*kB*Tb)) to avoid exp overflow)
            C_r = (Λ_2s_1s + Λ_α) / (Λ_2s_1s + Λ_α + β2)
            return C_r / H(co,x) * (β*(1-Xe) - nH(co,x)*α2*Xe^2)
        end
        co.Xe_Peebles_spline = _spline_integral(dXe_dx, x1, +20.0, Xe1, 1e-15)
    end

    return co.Xe_Peebles_spline(x) # TODO: spline the logarithm instead?
end

time_switch_Peebles(co::ΛCDM; Xe0::Real=0.999) = find_zero(x -> Xe_Saha_H_He(co, x) - Xe0, (-8.0, -7.0), rtol=1e-20, atol=1e-20)

time_reionization_H(co::ΛCDM)  = co.reionization ? -log(1 + co.z_reion_H)  : NaN
time_reionization_He(co::ΛCDM) = co.reionization ? -log(1 + co.z_reion_He) : NaN

function Xe_reionization(co::ΛCDM, x::Real)
    if !co.reionization
        return 0
    end

    y(z) = (1+z)^(3/2)
    dy_dz(z) = 3/2 * (1+z)^(1/2)
    Δy(z, Δz) = dy_dz(z) * Δz
    smoothstep(y, Δy, h) = h/2 * (1 + tanh(y / Δy))
    # TODO: replace heights with my hlines?
    Xe_reionization_total  = smoothstep(y(co.z_reion_H ) - y(z(x)), Δy(co.z_reion_H,  co.Δz_reion_H),  1+fHe(co.Yp))
    Xe_reionization_total += smoothstep(y(co.z_reion_He) - y(z(x)), Δy(co.z_reion_He, co.Δz_reion_He), 0+fHe(co.Yp))
    return Xe_reionization_total
end

# TODO: spline whole thing?
function Xe(co::ΛCDM, x::Real)
    if isnan(co.x_switch_Peebles)
        co.x_switch_Peebles = time_switch_Peebles(co)
    end
    Xe0 = Xe_Saha_H_He(co, co.x_switch_Peebles) # ≈ start from corresponding value from Saha

    if x < co.x_switch_Peebles
        Xe_total = Xe_Saha_H_He(co, x)
    else
        Xe_total = Xe_Peebles(co, x, co.x_switch_Peebles, Xe0)
    end

    Xe_total += Xe_reionization(co, x)

    return Xe_total
end

ne(co::ΛCDM, x::Real) = nH(co,x) * Xe(co,x)

dτ(co::ΛCDM, x::Real) = -ne(co,x) * σT * c / H(co,x)

function τ(co::ΛCDM, x::Real; derivative::Integer=0)
    if isnothing(co.τ_spline)
        co.τ_spline = _spline_integral((x, τ) -> dτ(co, x), 0.0, -20.0, 0.0, 1e-15)
    end
    return co.τ_spline(x, Val{derivative})
end

d2τ(co::ΛCDM, x::Real) = τ(co, x; derivative=2)
d3τ(co::ΛCDM, x::Real) = τ(co, x; derivative=3)

  g(co::ΛCDM, x::Real) = -dτ(co,x) * exp(-τ(co,x))
 dg(co::ΛCDM, x::Real) = -d2τ(co,x)*exp(-τ(co,x)) + dτ(co,x)^2*exp(-τ(co,x))
d2g(co::ΛCDM, x::Real) = (d3τ(co,x)*dτ(co,x)-d2τ(co,x)^2) / dτ(co,x)^2 * g(co,x) + d2τ(co,x)/dτ(co,x)*dg(co,x) - d2τ(co,x)*g(co,x) - dτ(co,x)*dg(co,x)

time_last_scattering_surface(co::ΛCDM) = find_zero(x -> τ(co, x) - 1, (-20, 0))
time_recombination(co::ΛCDM) = find_zero(x -> Xe(co, x) - 0.1, (-20, -3))

function sound_horizon(co::ΛCDM, x::Real)
    R(x) = 4*co.Ωγ0 / (3*co.Ωb0*a(x))
    cs(x) = c * √(R(x) / (3*(1+R(x))))
    ds_dx(x, s) = cs(x) / aH(co, x)
    x0 = -20.0
    s0 = cs(x0) / aH(co, x0)
    return _spline_integral(ds_dx, x0, +20.0, s0, 1e-10)(x)
end
