export Background
export H, aH, aH′, aH′′
export Ωγ, Ων, Ωb, Ωc, Ωk, ΩΛ, Ωr, Ωm, Ω
export t, t_radiation_matter_dominated, t_radiation_dominated
export η, η_radiation_matter_dominated, η_radiation_dominated
export χ, r, dL, dA
export x_equality_rm, x_equality_mΛ, x_acceleration

struct Background
    par::Parameters
    η::ODESolution # conformal time
    t::ODESolution # cosmic    time

    function Background(par::Parameters; x1::Float64=-20.0, x2::Float64=+20.0)
        η = integrate_η(par, x1, x2)
        t = integrate_t(par, x1, x2)
        new(par, η, t)
    end
end

aeq(par::Parameters) = a(x_equality_rm(par)) # for convenience below

η_radiation_matter_dominated(par::Parameters, x) = 2 / (par.H0*√(par.Ωm0)) * (√(a(x)+aeq(par)) - √(aeq(par))) # analytical solution with Ωk=ΩΛ=0
η_radiation_dominated(par::Parameters, x) = 1 / aH(par, x) # analytical solution with Ωm=Ωk=ΩΛ=0
η(bg::Background, x) = bg.η(x) # numerical solution

t_radiation_matter_dominated(par::Parameters, x) = 2 / (3*par.H0*√(par.Ωm0)) * (√(a(x)+aeq(par)) * (a(x)-2*aeq(par)) + 2*aeq(par)^(3/2)) # analytical solution with Ωk=ΩΛ=0
t_radiation_dominated(par::Parameters, x) = 1 / (2*H(par, x)) # analytical solution with Ωm=Ωk=ΩΛ=0
t(bg::Background, x) = bg.t(x) # numerical solution

# Integrate conformal time from x1 to x2
function integrate_η(par::Parameters, x1, x2)
    dη_dx(x, η) = 1 / aH(par, x)
    η1 = par.Ωm0 == 0 ? η_radiation_dominated(par, x1) : η_radiation_matter_dominated(par, x1)
    return solve(ODEProblem((η,_,x) -> dη_dx(x,η), η1, (x1, x2)), Tsit5(); abstol=1e-8, reltol=1e-8)
end

function integrate_t(par::Parameters, x1, x2)
    dt_dx(x, η) = 1 / H(par, x)
    t1 = par.Ωm0 == 0 ? t_radiation_dominated(par, x1) : t_radiation_matter_dominated(par, x1)
    return solve(ODEProblem((t,_,x) -> dt_dx(x,t), t1, (x1, x2)), Tsit5(); abstol=1e-8, reltol=1e-8)
end

# Friedmann equation
E(par::Parameters, x) = par.Ωr0/a(x)^4 + par.Ωm0/a(x)^3 + par.Ωk0/a(x)^2 + par.ΩΛ0 # H = H0*√(E)
H(par::Parameters, x) = par.H0 * √(E(par, x)) # cosmic    Hubble parameter
aH(par::Parameters, x) = a(x) * H(par, x)     # conformal Hubble parameter
aH′(par::Parameters, x) = ForwardDiff.derivative(x ->  aH(par, x), x)
aH′′(par::Parameters, x) = ForwardDiff.derivative(x -> aH′(par, x), x)

# density parameters (relative to critical density *at the time*)
# computed using Ωs = ρs/ρcrit = ρs/ρcrit0 * ρcrit0/ρcrit = Ωs0 * H0^2/H^2 = Ωs0 / E
Ωγ(par::Parameters, x) = par.Ωγ0 / a(x)^4 / E(par, x)
Ων(par::Parameters, x) = par.Ων0 / a(x)^4 / E(par, x)
Ωb(par::Parameters, x) = par.Ωb0 / a(x)^3 / E(par, x)
Ωc(par::Parameters, x) = par.Ωc0 / a(x)^3 / E(par, x)
Ωk(par::Parameters, x) = par.Ωk0 / a(x)^2 / E(par, x)
ΩΛ(par::Parameters, x) = par.ΩΛ0          / E(par, x)
Ωr(par::Parameters, x) = Ωγ(par, x) + Ων(par, x)
Ωm(par::Parameters, x) = Ωb(par, x) + Ωc(par, x)
Ω( par::Parameters, x) = Ωr(par, x) + Ωm(par, x) + Ωk(par, x) + ΩΛ(par, x)

# time of equality between different species and acceleration onset (as x = log(a))
x_equality_rm(par::Parameters) = x(par.Ωr0/par.Ωm0)
x_equality_mΛ(par::Parameters) = x(par.Ωm0/par.ΩΛ0) / 3
x_acceleration(par::Parameters) = find_zero(x -> aH′(par, x), (-20, +20))

# checks whether the Hubble parameter becomes zero (and complex) on the integration interval (x1, x2)
# EXAMPLE: has_turnaround(Parameters(Ωb0=0, Ωc0=0.2, Ωk0=-0.9)) == true
function has_turnaround(par::Parameters; x1=-20.0, x2=+20.0)
    # negative at the endpoints?
    if E(par, x1) < 0 || E(par, x2) < 0
        return true
    end

    # negative between the endpoints?
    # check the values at the stationary points of E(x), defined by roots of 2nd degree polynomial
    # dE/dx = -a^(-2) * [ 4*Ωr0*a^(-2) + 3*Ωm0*a^(-1) + 2*Ωk0 ] = 0
    ainv1, ainv2 = quadroots(4*par.Ωr0, 3*par.Ωm0, 2*par.Ωk0)
    if !isnan(ainv1)
        a1, a2 = 1/ainv1, 1/ainv2
        if (a1 >= 0 && E(par, x(a1)) < 0) || (a2 >= 0 && E(par, x(a2)) < 0)
            return true
        end
    end

    # never negative, so always positive
    return false
end

# conformal distance
χ(bg::Background, x) = c * (η(bg,0) - η(bg,x))

# radial coordinate (of light emitted at x)
#   note: ALL three expressions
#   1. r(Ωk0 = 0) = χ
#   2. r(Ωk0 < 0) = χ *  sin(√(-Ωk0)*H0*χ/c) / (√(-Ωk0)*H0*χ/c)
#   3. r(Ωk0 > 0) = χ * sinh(√(+Ωk0)*H0*χ/c) / (√(+Ωk0)*H0*χ/c)
#   can be written as the real-valued function (with intermediate complex results)
#      r(Ωk0)     = χ * sinc(√(-Ωk0)*H0*χ/c/π)
#   because sinc(x) = sin(π*x) / (π*x) -> 1 as x -> 0, and sinh(x) = -i * sin(i*x)
r(bg::Background, x) = χ(bg, x) * real(sinc(√(complex(-bg.par.Ωk0)) * bg.par.H0 * χ(bg, x) / c / π)) # in Julia, sinc(x) = sin(π*x) / (π*x), so divide argument by π!

# angular diameter distance and luminosity distance (of light emitted at x)
dA(bg::Background, x) = r(bg, x) * a(x)
dL(bg::Background, x) = r(bg, x) / a(x)
