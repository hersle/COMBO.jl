#struct Recombination{TXe <: ODESolution, Tτ <: ODESolution, Ts <: ODESolution}
struct Recombination
    bg::Background

    xswitch::Float64
    Xe_spline::ODESolution # free electron fraction
    τ_spline::ODESolution # optical depth
    sound_horizon_spline::ODESolution # sound horizon

    #Recombination(bg::Background, xswitch::Float64) where {Tη <: ODESolution, Tt <: ODESolution} = new{Tη,Tt}(par, η_spline, t_spline)
end

function Recombination(bg::Background; x0=-20.0)
    xswitch = time_switch_Peebles(bg.par) # TODO: why does this allocate?
    Xe_spline = spline_Xe(bg.par, xswitch)

    Xe(x) = (x < xswitch ? Xe_Saha_H_He(bg.par, x) : Xe_spline(x)) + Xe_reionization(bg.par, x)
    ne(x) = nH(bg.par,x) * Xe(x)
    dτ_dx(x) = -ne(x) * σT * c / H(bg.par,x)
    τ_spline = solve(ODEProblem((τ,_,x) -> dτ_dx(x), 0.0, (0.0, x0)), Tsit5(); abstol=1e-8, reltol=1e-8)

    R(x) = 4*bg.par.Ωγ0 / (3*bg.par.Ωb0*a(x))
    cs(x) = c * √(R(x) / (3*(1+R(x))))
    ds_dx(x, s) = cs(x) / aH(bg.par, x)
    s0 = cs(x0) / aH(bg.par, x0)
    sound_horizon_spline = solve(ODEProblem((s,_,x) -> ds_dx(x,s), s0, (x0, 0.0)), Tsit5(); abstol=1e-8, reltol=1e-8) # TODO: rename s

    Recombination(bg, xswitch, Xe_spline, τ_spline, sound_horizon_spline)
end

Base.broadcastable(rec::Recombination) = Ref(rec)

ρcrit(par::Parameters, x::Real) = 3 * H(par, x)^2 / (8 * π * G)
ρb(par::Parameters, x::Real) = Ωb(par, x) * ρcrit(par, x)
nb(par::Parameters, x::Real) = ρb(par, x) / mH
nH(par::Parameters, x::Real) = (1-par.Yp) * nb(par, x)
nHe(par::Parameters, x::Real) = par.Yp/4 * nb(par, x)
Tγ(par::Parameters, x::Real) = par.Tγ0 / a(x)
fHe(Yp::Real) = Yp / (4*(1-Yp))

function Xe_Saha_H(par::Parameters, x::Real)
    Tb = Tγ(par, x)
    λe = h / √(2*π*me*kB*Tb)
    a = 1
    b = 1 / λe^3 * exp(-EH1ion/(kB*Tb)) / nH(par,x)
    c = -b

    # when b >> 1, the quadratic equation solution is
    #   (-b + √(b^2+4*b)) / 2
    # = b/2 * (-1 + √(1+4/b))
    # ≈ b/2 * (-1 + 1 + 2/b - 2/b^2)  (b >> 1)
    # = 1 - 1/b
    return b < 1e10 ? quadroots(a, b, c)[2] : 1 - 1/b # choose Taylor expansion for large b
end

function Xe_Saha_H_He(par::Parameters, x::Real; tol::Float64=1e-15, maxiters::Int=10000)
    # Fixed-point iterate for Xe, starting with the fast H-only Saha equation as the initial guess
    # If Xe0 is small, the iterations converges very slowly, but then there is almost no He anyway, so we just take the H-only solution
    Xe0 = Xe_Saha_H(par, x)

    # To find Xe = XH+ + Yp/(4*(1-Yp)) * (XHe+ + 2*XHe++),
    # begin with an initial guess for Xe and iteratively solve the system of Saha equations
    # (1) ne * XHe+ / (1 - XHe+ - XHe++) = 2 / λe^3 * exp(-EHe1ion / (kB*Tb))
    # (2) ne * XHe++ / XHe+ = 4 / λe^3 * exp(-EHe2ion/(kB*Tb))
    # (3) ne * XH+ / (1-XH+) = 1 / λe^3 * exp(-EH1ion/(kB*Tb))
    # for {XH+, XHe+, XHe++} = {XH1, XHe1, XHe2}.
    function Xe_Saha_H_He_fixed_point(Xe)
        Tb = Tγ(par, x)
        λe = h / √(2*π*me*kB*Tb)
        ne = Xe * nH(par, x)
        R1 = 1 / λe^3 * exp(-EH1ion /(kB*Tb)) / ne
        R2 = 2 / λe^3 * exp(-EHe1ion/(kB*Tb)) / ne
        R3 = 4 / λe^3 * exp(-EHe2ion/(kB*Tb)) / ne
        XH1  = 1 / (1 + 1/R1)
        XHe1 = 1 / (1 + 1/R2 + R3)
        XHe2 = 1 / (1 + 1/R3 + 1/(R2*R3))
        return XH1 + fHe(par.Yp) * (XHe1 + 2*XHe2)
    end

    return Xe0 < 0.5 ? Xe0 : fixed_point_iterate(Xe_Saha_H_He_fixed_point, Xe0; tol=tol)
end

# @code_warntype on Xe(co,x) says that this can return Any, unless its return type is explicitly stated (due to Union{...., Nothing}?)
# TODO: handle in a better way?
function Xe_Peebles_spline(par::Parameters, x1::Float64, Xe1::Float64)
    function dXe_dx(x, Xe)
        Tb = Tγ(par,x) # K # TODO: assumptions? separate baryon evolution?
        n_1s = (1-Xe) * nH(par,x) # 1/m^3
        Λ_2s_1s = 8.227 # 1/s
        Λ_α = H(par,x) * (3*EH1ion/(ħ*c))^3 / ((8*π)^2 * n_1s) # 1/s
        ϕ2 = 0.448 * log(EH1ion/(kB*Tb)) # dimensionless
        α2 = 64*π / √(27*π) * (α/me)^2 * √(EH1ion/(kB*Tb)) * ϕ2 * ħ^2/c # m^3/s
        λe = √(h^2 / (2*π*me*kB*Tb)) # thermal de Broglie wavelength
        β  = α2 / λe^3 * exp(-EH1ion/(kB*Tb))
        β2 = α2 / λe^3 * exp(-EH1ion/(4*kB*Tb)) # 1/s (compute this instead of β2 = β * exp(3*EH1ion/(4*kB*Tb)) to avoid exp overflow)
        C_r = (Λ_2s_1s + Λ_α) / (Λ_2s_1s + Λ_α + β2)
        return C_r / H(par,x) * (β*(1-Xe) - nH(par,x)*α2*Xe^2)
    end
    #return _spline_integral(dXe_dx, x1, +20.0, Xe1; name="free electron fraction Xe")
    return solve(ODEProblem((Xe,_,x) -> dXe_dx(x, Xe), Xe1, (x1, 0.0)), Tsit5(); abstol=1e-8, reltol=1e-8)
end

function time_switch_Peebles(par::Parameters)::Float64
    return find_zero(x -> Xe_Saha_H_He(par, x) - 0.999, (-20.0, +20.0), rtol=1e-20, atol=1e-20)
end

time_reionization_H(par::Parameters)  = par.reionization ? -log(1 + par.z_reion_H)  : NaN
time_reionization_He(par::Parameters) = par.reionization ? -log(1 + par.z_reion_He) : NaN

function Xe_reionization(par::Parameters, x::Real)
    if !par.reionization
        return 0.0
    end

    y(z) = (1+z)^(3/2)
    dy_dz(z) = 3/2 * (1+z)^(1/2)
    Δy(z, Δz) = dy_dz(z) * Δz
    smoothstep(y, Δy, h) = h/2 * (1 + tanh(y / Δy))
    Xe_reionization_total  = smoothstep(y(par.z_reion_H ) - y(z(x)), Δy(par.z_reion_H,  par.Δz_reion_H),  1+fHe(par.Yp))
    Xe_reionization_total += smoothstep(y(par.z_reion_He) - y(z(x)), Δy(par.z_reion_He, par.Δz_reion_He), 0+fHe(par.Yp))
    return Xe_reionization_total
end

# TODO: check type stability here, now that I do branching on Saha/Peebles
function spline_Xe(par::Parameters, xswitch::Real; x1::Float64=-20.0)
    return Xe_Peebles_spline(par, xswitch, Xe_Saha_H_He(par, xswitch)) # start Peebles from last value of Saha
end

# TODO: don't duplicate in constructor
Xe(rec::Recombination, x) = (x < rec.xswitch ? Xe_Saha_H_He(rec.bg.par, x) : rec.Xe_spline(x)) + Xe_reionization(rec.bg.par, x)

  τ(rec::Recombination, x) = rec.τ_spline(x)
 dτ(rec::Recombination, x) = ForwardDiff.derivative(x ->  τ(rec, x), x)
d2τ(rec::Recombination, x) = ForwardDiff.derivative(x -> dτ(rec, x), x)

  g(rec::Recombination, x) = -dτ(rec, x) * exp.(-τ(rec, x))
 dg(rec::Recombination, x) = ForwardDiff.derivative(x ->  g(rec, x), x)
d2g(rec::Recombination, x) = ForwardDiff.derivative(x -> dg(rec, x), x)

sound_horizon(rec::Recombination, x) = rec.sound_horizon_spline(x)

time_decoupling(rec::Recombination) = find_zero(x -> dτ(rec,x)^2 - d2τ(rec,x) - 0.0, (-20.0, -3.0)) # equivalent to dg=0 without the exponential; exclude reionization for x > -3
time_recombination(rec::Recombination) = find_zero(x -> Xe(rec,x) - 0.1, (-20.0, -3.0)) # exclude reionization for x > -3
