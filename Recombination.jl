ρcrit(co::ΛCDM, x::Real) = 3 * H(co, x)^2 / (8 * π * G)
ρb(co::ΛCDM, x::Real) = Ωb(co, x) * ρcrit(co, x)
nb(co::ΛCDM, x::Real) = ρb(co, x) / mH
nH(co::ΛCDM, x::Real) = (1-co.Yp) * nb(co, x)
nHe(co::ΛCDM, x::Real) = co.Yp/4 * nb(co, x)
Tγ(co::ΛCDM, x::Real) = co.Tγ0 / a(x)
fHe(Yp::Real) = Yp / (4*(1-Yp))

function Xe_Saha_H(co::ΛCDM, x::Real)
    Tb = Tγ(co, x)
    λe = h / √(2*π*me*kB*Tb)
    a = 1
    b = 1 / λe^3 * exp(-EH1ion/(kB*Tb)) / nH(co,x)
    c = -b

    # when b >> 1, the quadratic equation solution is
    #   (-b + √(b^2+4*b)) / 2
    # = b/2 * (-1 + √(1+4/b))
    # ≈ b/2 * (-1 + 1 + 2/b - 2/b^2)  (b >> 1)
    # = 1 - 1/b
    return b < 1e10 ? quadroots(a, b, c)[2] : 1 - 1/b # choose Taylor expansion for large b
end

# To find Xe = XH+ + Yp/(4*(1-Yp)) * (XHe+ + 2*XHe++),
# begin with an initial guess for Xe and iteratively solve the system of Saha equations
# (1) ne * XHe+ / (1 - XHe+ - XHe++) = 2 / λe^3 * exp(-EHe1ion / (kB*Tb))
# (2) ne * XHe++ / XHe+ = 4 / λe^3 * exp(-EHe2ion/(kB*Tb))
# (3) ne * XH+ / (1-XH+) = 1 / λe^3 * exp(-EH1ion/(kB*Tb))
# for {XH+, XHe+, XHe++} = {XH1, XHe1, XHe2}.
function Xe_Saha_H_He_fixed_point(co, x, Xe)
    Tb = Tγ(co, x)
    λe = h / √(2*π*me*kB*Tb)
    ne = Xe * nH(co, x)
    R1 = 1 / λe^3 * exp(-EH1ion /(kB*Tb)) / ne
    R2 = 2 / λe^3 * exp(-EHe1ion/(kB*Tb)) / ne
    R3 = 4 / λe^3 * exp(-EHe2ion/(kB*Tb)) / ne
    XH1  = 1 / (1 + 1/R1)
    XHe1 = 1 / (1 + 1/R2 + R3)
    XHe2 = 1 / (1 + 1/R3 + 1/(R2*R3))
    return XH1 + fHe(co.Yp) * (XHe1 + 2*XHe2)
end

function Xe_Saha_H_He(co::ΛCDM, x::Real; tol::Float64=1e-15, maxiters::Int=10000)
    # Fixed-point iterate for Xe, starting with the fast H-only Saha equation as the initial guess
    # If Xe0 is small, the iterations converges very slowly, but then there is almost no He anyway, so we just take the H-only solution
    Xe0 = Xe_Saha_H(co, x)
    return Xe0 < 0.5 ? Xe0 : fixed_point_iterate(Xe -> Xe_Saha_H_He_fixed_point(co, x, Xe), Xe0; tol=tol)
end

# @code_warntype on Xe(co,x) says that this can return Any, unless its return type is explicitly stated (due to Union{...., Nothing}?)
# TODO: handle in a better way?
function Xe_Peebles_spline(co::ΛCDM, x1::Float64, Xe1::Float64)
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
    return _spline_integral(dXe_dx, x1, +20.0, Xe1; name="free electron fraction Xe")
end

function time_switch_Peebles(co::ΛCDM)::Float64
    return find_zero(x -> Xe_Saha_H_He(co, x) - 0.999, (-20.0, +20.0), rtol=1e-20, atol=1e-20)
end

time_reionization_H(co::ΛCDM)  = co.reionization ? -log(1 + co.z_reion_H)  : NaN
time_reionization_He(co::ΛCDM) = co.reionization ? -log(1 + co.z_reion_He) : NaN

function Xe_reionization(co::ΛCDM, x::Real)
    if !co.reionization
        return 0.0
    end

    y(z) = (1+z)^(3/2)
    dy_dz(z) = 3/2 * (1+z)^(1/2)
    Δy(z, Δz) = dy_dz(z) * Δz
    smoothstep(y, Δy, h) = h/2 * (1 + tanh(y / Δy))
    Xe_reionization_total  = smoothstep(y(co.z_reion_H ) - y(z(x)), Δy(co.z_reion_H,  co.Δz_reion_H),  1+fHe(co.Yp))
    Xe_reionization_total += smoothstep(y(co.z_reion_He) - y(z(x)), Δy(co.z_reion_He, co.Δz_reion_He), 0+fHe(co.Yp))
    return Xe_reionization_total
end

function Xe(co::ΛCDM, x::Real; x1::Float64=-20.0)
    if isnothing(co.Xe_spline)
        xswitch = time_switch_Peebles(co) # TODO: why does this allocate?

        # Peebles equation points
        x2, Xe2spl = Xe_Peebles_spline(co, xswitch, Xe_Saha_H_He(co, xswitch)) # start Peebles from last value of Saha

        # Saha equation points
        x1 = range(x1, xswitch, length=length(x2)) # use as many points as Peebles; don't duplicate xswitch
        Xe1spl = spline(x1, Xe_Saha_H_He.(co, x1))

        # merge Saha and Peebles
        x12, co.Xe_spline = splinejoin(x1, x2, Xe1spl, Xe2spl) # spline is WITHOUT reionization (the spline points do not resolve reionization)!

        # Reionization points
        if co.reionization
            dx(z) = -1/(1+z)
            xH  = time_reionization_H(co)
            xHe = time_reionization_He(co)
            dxH  = 1/(1+co.z_reion_H) # = |x′(z)|
            dxHe = 1/(1+co.z_reion_He)
            x3H  = range(xH-5*dxH, xH+5*dxH, length=100) # resolve one reionization with 100 points extending ±5dx from central x (in addition to Saha/Peebles points in x12)
            x3He = range(xHe-5*dxHe, xHe+5*dxHe, length=100)
            x3 = unique(sort(vcat(x3H, x3He)))
        else
            x3 = []
        end

        # merge (Saha and Peebles) and reionization
        x123 = unique(sort(vcat(x12, x3)))
        co.Xe_spline = spline(x123, co.Xe_spline(x123) .+ Xe_reionization.(co, x123))
    end

    return co.Xe_spline(x)
end

ne(co::ΛCDM, x::Real) = nH(co,x) * Xe(co,x)

dτ(co::ΛCDM, x::Real) = -ne(co,x) * σT * c / H(co,x) # TODO: spline!!

function τ(co::ΛCDM, x::Real; deriv::Integer=0)
    if isnothing(co.τ_spline)
        _, co.τ_spline = _spline_integral((x, τ) -> dτ(co, x), 0.0, -20.0, 0.0; name="optical depth τ")
    end
    return deriv == 0 ? co.τ_spline(x) : derivative(co.τ_spline, x; nu=deriv)
end

d2τ(co::ΛCDM, x::Real) = τ(co, x; deriv=2)
d3τ(co::ΛCDM, x::Real) = τ(co, x; deriv=3)

function g(co::ΛCDM, x::Real; deriv::Integer=0)
    if isnothing(co.g_spline)
        τ(co, 0.0) # trigger spline computation
        xs = extendx(splinex(co.τ_spline), 3)
        co.g_spline = spline(xs, @. -dτ(co,xs) * exp(-τ(co,xs)))
    end
    return deriv == 0 ? co.g_spline(x) : derivative(co.g_spline, x; nu=deriv)
end

 dg(co::ΛCDM, x::Real) = g(co, x; deriv=1)
d2g(co::ΛCDM, x::Real) = g(co, x; deriv=2)

time_decoupling(co::ΛCDM) = find_zero(x -> dτ(co,x)^2 - d2τ(co,x) - 0.0, (-20.0, -3.0)) # equivalent to dg=0 without the exponential; exclude reionization for x > -3
time_recombination(co::ΛCDM) = find_zero(x -> Xe(co, x) - 0.1, (-20.0, -3.0)) # exclude reionization for x > -3

function sound_horizon(co::ΛCDM, x::Real)
    if isnothing(co.sound_horizon_spline)
        R(x) = 4*co.Ωγ0 / (3*co.Ωb0*a(x))
        cs(x) = c * √(R(x) / (3*(1+R(x))))
        ds_dx(x, s) = cs(x) / aH(co, x)
        x0 = -20.0
        s0 = cs(x0) / aH(co, x0)
        _, co.sound_horizon_spline = _spline_integral(ds_dx, x0, +20.0, s0; name="sound horizon s")
    end
    return co.sound_horizon_spline(x)
end
