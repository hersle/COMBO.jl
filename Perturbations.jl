using OffsetArrays # arbitrary-indexed vectors

export PerturbationMode, PerturbationModes
export δc, δb, vc, vb, Φ, Ψ, Θl, Nl, ΘPl, ST, SE
export x_horizon_entry

const polarization = true
const neutrinos = true
const lmax = 10 # use >= 6 without neutrinos/polarization; >= 10 with neutrinos/polarization
@assert lmax >= 4 # equations are ambiguous with lmax <= 3 (how to handle l<2, l=2 and l=lmax?)

# variable index map (1-based),
const i_δc = 1
const i_vc = 2
const i_δb = 3
const i_vb = 4
const i_Φ  = 5
const i_Nl(l) = 6 + l
const i_Θl(l) = i_Nl(lmax) + 1 + l
const i_ΘPl(l) = i_Θl(lmax) + 1 + l
const i_max = i_ΘPl(lmax) # last variable of system

struct PerturbationMode
    rec::Recombination
    k::Float64 # wavenumber
    y::ODESolution # vector solution of entire system

    function PerturbationMode(rec::Recombination, k::Float64; kwargs...)
        new(rec, k, integrate_perturbation_mode(rec, k; kwargs...))
    end
end

function PerturbationModes(rec::Recombination, ks::AbstractVector; kwargs...)
    perturbs = Vector{PerturbationMode}(undef, length(ks))
    Threads.@threads for i in 1:length(ks)
        perturbs[i] = PerturbationMode(rec, ks[i]; kwargs...)
    end
    return perturbs
end

function perturbations_initial_conditions(par::Parameters, η::ODESolution, τ::ODESolution, x0, k)
    fν = par.Ων0 / par.Ωr0
    Ψ  = -1 / (3/2 + 2*fν/5)
    τ′ = ForwardDiff.derivative(τ, x0)

    y = Vector{Float64}(undef, i_max)
    y[i_δc] = y[i_δb] = -3/2 * Ψ
    y[i_vc] = y[i_vb] = -k*c/(2*aH(par,x0)) * Ψ
    y[i_Φ] = -(1 + 2*fν/5) * Ψ
    y[i_Θl(0)] = -1/2 * Ψ
    y[i_Θl(1)] = -c*k / (3*aH(par,x0)) * y[i_Θl(0)]
    y[i_Θl(2)] = (polarization ? -8/15 : -20/45) * c*k / (aH(par,x0)*τ′) * y[i_Θl(1)]
    for l in 3:lmax
        y[i_Θl(l)] = -l/(2*l+1) * c*k/(aH(par,x0)*τ′) * y[i_Θl(l-1)] # recursive relation
    end

    if polarization
        y[i_ΘPl(0)] = 5/4 * y[i_Θl(2)]
        y[i_ΘPl(1)] = -c*k/(4*aH(par,x0)*τ′) * y[i_Θl(2)]
        y[i_ΘPl(2)] = 1/4 * y[i_Θl(2)]
        for l in 3:lmax
            y[i_ΘPl(l)] = -l/(2*l+1) * c*k/(aH(par,x0)*τ′) * y[i_ΘPl(l-1)] # recursive relation
        end
    else
        y[i_ΘPl(0):i_ΘPl(lmax)] .= 0.0 # pin to zero
    end

    if neutrinos
        y[i_Nl(0)] = -1/2 * Ψ
        y[i_Nl(1)] = +c*k/(6*aH(par,x0)) * Ψ
        y[i_Nl(2)] = +(c*k*a(x0)/par.H0)^2 / (30*par.Ωr0) * Ψ # expanded to avoid 1/Ων0 with Ων0=0
        for l in 3:lmax
            y[i_Nl(l)] = c*k/((2*l+1)*aH(par,x0)) * y[i_Nl(l-1)]
        end
    else
        y[i_Nl(0):i_Nl(lmax)] .= 0.0 # pin to zero
    end

    return y
end

function integrate_perturbation_mode(rec::Recombination, k; kwargs...)
    return integrate_perturbation_mode(rec.bg.par, rec.bg.η, rec.τ, k; kwargs...)
end

# putting callables for η and τ behind a function barrier removes allocations!
function integrate_perturbation_mode(par::Parameters, η::ODESolution, τ::ODESolution, k; x1=-20.0, x2=0.0, verbose=true)
    function dy_dx!(y′, y, _, x)
        # 1) unpack variables from vector
        δc = y[i_δc]
        vc = y[i_vc]
        δb = y[i_δb]
        vb = y[i_vb]
        Φ  = y[i_Φ]
        Nl = OffsetVector(view(y, i_Nl(0):i_Nl(lmax)), 0:lmax) # shift indexing from Julia's default 1:lmax+1 to more natural 0:lmax
        Θl = OffsetVector(view(y, i_Θl(0):i_Θl(lmax)), 0:lmax)
        ΘPl = OffsetVector(view(y, i_ΘPl(0):i_ΘPl(lmax)), 0:lmax)

        # 2) pre-compute commonly combined auxilliary quantities
        ck_aH = c*k / aH(par,x)
        τ′ = ForwardDiff.derivative(τ, x)
        R = 4*par.Ωγ0 / (3*par.Ωb0*a(x))
        Ψ = -Φ - 12 * (par.H0 / (c*k*a(x)))^2 * (par.Ωγ0*Θl[2]+par.Ων0*Nl[2])
        Π = Θl[2] + ΘPl[0] + ΘPl[2]
        Φ′ = Ψ - ck_aH^2/3*Φ + (par.H0/aH(par,x))^2/2 * (par.Ωc0/a(x)*δc + par.Ωb0/a(x)*δb + 4*par.Ωγ0/a(x)^2*Θl[0] + 4*par.Ων0/a(x)^2*Nl[0])

        # 3) compute derivatives while packing them into a vector
        y′[i_δc]    = ck_aH*vc - 3*Φ′
        y′[i_δb]    = ck_aH*vb - 3*Φ′
        y′[i_vc]    = -vc - ck_aH*Ψ
        y′[i_vb]    = -vb - ck_aH*Ψ + τ′*R*(3*Θl[1]+vb)
        y′[i_Φ]     = Φ′

        y′[i_Θl(0)] = -ck_aH*Θl[1] - Φ′
        y′[i_Θl(1)] =  ck_aH/3 * (Θl[0]-2*Θl[2]+Ψ) + τ′*(Θl[1]+vb/3)
        for l in 2:lmax-1
            y′[i_Θl(l)] = ck_aH/(2*l+1) * (l*Θl[l-1] - (l+1)*Θl[l+1]) + τ′*(Θl[l]-Π/10*δ(l,2))
        end
        y′[i_Θl(lmax)] = ck_aH*Θl[lmax-1] - (lmax+1)/(aH(par,x)*η(x))*Θl[lmax] + τ′*Θl[lmax] # 2nd term: their η is my c*η

        if polarization
            y′[i_ΘPl(0)] = -ck_aH*ΘPl[1] + τ′*(ΘPl[0]-Π/2)
            for l in 1:lmax-1
                y′[i_ΘPl(l)] = ck_aH/(2*l+1) * (l*ΘPl[l-1] - (l+1)*ΘPl[l+1]) + τ′*(ΘPl[l]-Π/10*δ(l,2))
            end
            y′[i_ΘPl(lmax)] = ck_aH*ΘPl[lmax-1] - (lmax+1)/(aH(par,x)*η(x))*ΘPl[lmax] + τ′*ΘPl[lmax]
        else
            y′[i_ΘPl(0):i_ΘPl(lmax)] .= 0.0 # pin to zero
        end

        if neutrinos
            y′[i_Nl(0)] = -ck_aH*Nl[1] - Φ′
            y′[i_Nl(1)] =  ck_aH/3 * (Nl[0] - 2*Nl[2] + Ψ)
            for l in 2:lmax-1
                y′[i_Nl(l)] = ck_aH/(2*l+1) * (l*Nl[l-1] - (l+1)*Nl[l+1])
            end
            y′[i_Nl(lmax)] = ck_aH*Nl[lmax-1] - (lmax+1)/(aH(par,x)*η(x)) * Nl[lmax]
        else
            y′[i_Nl(0):i_Nl(lmax)] .= 0.0 # pin to zero
        end

        return nothing # integration uses y′ in-place
    end


    time = @elapsed begin
    y1 = perturbations_initial_conditions(par, η, τ, x1, k)
    prob = ODEProblem{true}(dy_dx!, y1, (x1, x2))
    sol = solve(prob, KenCarp4(autodiff=false); abstol=1e-5, reltol=1e-5)
    end
    if verbose
        println("Integrated perturbation mode k=$(k*Mpc)/Mpc with $(length(sol.t)) points in $(time) seconds")
    end
    return sol
end

# raw quantities (from integration)
Φ(perts::PerturbationMode, x) = perts.y(x)[i_Φ]
δc(perts::PerturbationMode, x) = perts.y(x)[i_δc]
δb(perts::PerturbationMode, x) = perts.y(x)[i_δb]
vc(perts::PerturbationMode, x) = perts.y(x)[i_vc]
vb(perts::PerturbationMode, x) = perts.y(x)[i_vb]
Θl(perts::PerturbationMode, x, l) = perts.y(x)[i_Θl(l)]
Nl(perts::PerturbationMode, x, l) = perts.y(x)[i_Nl(l)]
ΘPl(perts::PerturbationMode, x, l) = perts.y(x)[i_ΘPl(l)]

# "composite" quantities (form raw quantities)
function Ψ(perts::PerturbationMode, x; deriv=0)
    par = perts.rec.bg.par
    k = perts.k
    return -Φ(perts, x) - 12*par.H0^2/(c*k*a(x))^2 * (par.Ωγ0*Θl(perts, x, 2) + par.Ων0*Nl(perts, x, 2))
end
Π(perts::PerturbationMode, x) = Θl(perts, x, 2) + ΘPl(perts, x, 0) + ΘPl(perts, x, 2)

# source functions (for milestone 4)
S_SW(perts::PerturbationMode, x) = g(perts.rec,x) * (Θl(perts, x, 0) + Ψ(perts, x) + Π(perts,x)/4)
S_ISW(perts::PerturbationMode, x) = exp(-τ(perts.rec,x)) * ForwardDiff.derivative(x -> Ψ(perts,x) - Φ(perts,x), x)
S_Doppler(perts::PerturbationMode, x) = -1/(c*perts.k) * ForwardDiff.derivative(x -> aH(perts.rec.bg.par,x) * g(perts.rec,x) * vb(perts, x), x)
S_polarization(perts::PerturbationMode, x) = 3/(4*c^2*perts.k^2) * ForwardDiff.derivative(x -> aH(perts.rec.bg.par,x) * ForwardDiff.derivative(x -> aH(perts.rec.bg.par,x) * g(perts.rec,x) * Π(perts,x), x), x)
ST(perts::PerturbationMode, x) = S_SW(perts, x) + S_ISW(perts, x) + S_Doppler(perts, x) + S_polarization(perts, x)
SE(perts::PerturbationMode, x) = x == 0.0 ? 0.0 : 3 * g(perts.rec,x) * Π(perts,x) / (2*c*perts.k*(η(perts.rec.bg,0.0)-η(perts.rec.bg,x)))^2 # force SE(x=0) = 0 (otherwise gives Inf)

# horizon entry time
x_horizon_entry(bg::Background, k) = find_zero(x -> k*c*η(bg,x) - 1, (-20.0, 0.0))
