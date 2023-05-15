struct Background
    par::Parameters # TODO: pass by reference?
    η::ODESolution # conformal time
    t::ODESolution # cosmic    time

    function Background(par::Parameters; x1::Float64=-20.0, x2::Float64=+20.0)
        η = integrate_η(par, x1, x2)
        t = integrate_t(par, x1, x2)
        new(par, η, t)
    end
end

# Fix bg in vectorized calls like H.(bg, [-1.0, -2.0])
Base.broadcastable(bg::Background) = Ref(bg)

# Integrate conformal time from x1 to x2
function integrate_η(par, x1, x2)
    dη_dx(x, η) = 1 / aH(par, x)
    if par.Ωm0 > 0
        η1 = 2 / (par.H0*√(par.Ωm0)) * (√(a(x1)+aeq(par)) - √(aeq(par))) # analytical solution with Ωk=ΩΛ=0
    else
        η1 = 1 / aH(par, x1) # analytical solution with Ωm=Ωk=ΩΛ=0
    end
    η_spline = solve(ODEProblem((η,_,x) -> dη_dx(x,η), η1, (x1, x2)), Tsit5(); abstol=1e-8, reltol=1e-8)
end

function integrate_t(par, x1, x2)
    dt_dx(x, η) = 1 / H(par, x)
    if par.Ωm0 > 0
        t1 = 2 / (3*par.H0*√(par.Ωm0)) * (√(a(x1)+aeq(par)) * (a(x1)-2*aeq(par)) + 2*aeq(par)^(3/2)) # analytical solution with Ωk=ΩΛ=0
    else
        t1 = 1 / (2*H(par, x1)) # analytical solution with Ωm=Ωk=ΩΛ=0
    end
    t_spline = solve(ODEProblem((t,_,x) -> dt_dx(x,t), t1, (x1, x2)), Tsit5(); abstol=1e-8, reltol=1e-8)
end

η(bg::Background, x) = bg.η(x)
t(bg::Background, x) = bg.t(x)

# conformal distance
χ(bg::Background, x) = c * (η(bg,0) - η(bg,x))

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
r(bg::Background, x) = χ(bg, x) * real(sinc(√(complex(-bg.par.Ωk0)) * bg.par.H0 * χ(bg, x) / c / π)) # in Julia, sinc(x) = sin(π*x) / (π*x), so divide argument by π!

# angular diameter distance and luminosity distance (of light emitted at x)
dA(bg::Background, x) = r(bg, x) * a(x)
dL(bg::Background, x) = r(bg, x) / a(x)
