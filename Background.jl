struct Background
    par::Parameters # TODO: pass by reference?

    η_spline::ODESolution # conformal time
    t_spline::ODESolution # cosmic    time

    function Background(par::Parameters)
        # integrate η (TODO: build from Parameters struct?)
        dη_dx(x, η) = 1 / aH(par, x)
        x1, x2 = -20.0, +20.0
        aeq = par.Ωr0 / par.Ωm0
        if par.Ωm0 > 0
            η1 = 2 / (par.H0*√(par.Ωm0)) * (√(a(x1)+aeq) - √(aeq)) # anal expr with Ωk=ΩΛ=0
        else
            η1 = 1 / aH(par, x1) # anal expr with Ωm=Ωk=ΩΛ=0
        end
        _, η_spline = _spline_integral(dη_dx, x1, x2, η1; name="conformal time η")

        # integrate t
        dt_dx(x, η) = 1 / H(par, x)
        x1, x2 = -20.0, +20.0
        aeq = par.Ωr0 / par.Ωm0
        if par.Ωm0 > 0
            t1 = 2 / (3*par.H0*√(par.Ωm0)) * (√(a(x1)+aeq) * (a(x1)-2*aeq) + 2*aeq^(3/2)) # anal expr with Ωk=ΩΛ=0
        else
            t1 = 1 / (2*H(par, x1)) # anal expr with Ωm=Ωk=ΩΛ=0
        end
        _, t_spline = _spline_integral(dt_dx, x1, x2, t1; name="cosmic time t")

        new(par, η_spline, t_spline)
    end
end

# in vectorized calls, like H.(co, [1.0, 2.0]),
# broadcast the same cosmology to all scalar calls
# TODO: apply to recombination etc., too
Base.broadcastable(bg::Background) = Ref(bg)

η(bg::Background, x::Real) = bg.η_spline(x)
t(bg::Background, x::Real) = bg.t_spline(x)

# conformal distance
χ(bg::Background, x::Real) = c * (η(bg, 0) - η(bg, x))

# radial coordinate (of light emitted at x)
#   note: ALL three expressions
#   1. r(Ωk0 = 0) = χ
#   2. r(Ωk0 < 0) = χ *  sin(√(-Ωk0)*H0*χ/c) / (√(-Ωk0)*H0*χ/c)
#   3. r(Ωk0 > 0) = χ * sinh(√(+Ωk0)*H0*χ/c) / (√(+Ωk0)*H0*χ/c)
#   can be written as the real-valued function
#      r(Ωk0)     = χ * sinc(√(-Ωk0)*H0*χ/c/π)
#   using complex numbers,
#   because sinc(x) = sin(π*x) / (π*x) -> 1 as x -> 0,
#   and sinh(x) = -i * sin(i*x)
r(bg::Background, x::Real) = χ(bg, x) * real(sinc(√(complex(-bg.par.Ωk0)) * bg.par.H0 * χ(bg, x) / c / π)) # in Julia, sinc(x) = sin(π*x) / (π*x), so divide argument by π!

# angular diameter distance and luminosity distance (of light emitted at x)
dA(bg::Background, x::Real) = r(bg, x) * a(x)
dL(bg::Background, x::Real) = r(bg, x) / a(x)
