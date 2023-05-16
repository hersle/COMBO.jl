P_primordial(par::Parameters, k::Real) = 2*π^2/k^3 * par.As * (k/par.k_pivot)^(par.ns-1)
Δ(perturb::PerturbationMode, x::Real, k::Real) = c^2*k^2*Φ(perturb,x) / (3/2*perturb.rec.bg.par.Ωm0/a(x)*perturb.rec.bg.par.H0^2)
P(perturb::PerturbationMode, x::Real, k::Real) = abs(Δ(perturb,x,k))^2 * P_primordial(perturb.rec.bg.par, k)

# TODO:
# From Bolt.jl creator: Levin method for oscillatory Bessel integrals
# https://discourse.julialang.org/t/rfc-ann-oscillatoryintegralsode-jl-levin-method-ordinarydiffeq/55601

function integrate_trapz(f, x1, x2; step=nothing, length=nothing)
    @assert !isnothing(step) || !isnothing(length) "step or length must be specified"
    xs = range(x1, x2; step=step, length=length)
    return trapz(xs, f.(xs))
end

function integrate_adaptive(f, x1, x2; atol=0, rtol=1e-3, order=8) # TODO: rtol=1e-3 and order=8 gives good results!
    return quadgk(f, x1, x2; atol=atol, rtol=rtol, order=order)[1] # discard error
end

dΘl0_dx(l, k, x, S, η) = S(x, k) * jl(l, c*k*(η(0)-η(x)))
Θl0(l, k, S, η)::Float64 = integrate_trapz(x -> dΘl0_dx(l,k,x,S,η), -20.0, 0.0; step=0.02)

dΘEl0_dx(l::Integer, k, x, SE, η) = √((l+2)*(l+1)*(l+0)*(l-1)) * SE(x, k) * jl(l, c*k*(η(0.0)-η(x)))
ΘEl0(l::Integer, k, SE, η)::Float64 = integrate_trapz(x -> dΘEl0_dx(l,k,x,SE,η), -20.0, 0.0; step=0.02)

dCl_dk_generic(l,k,ΘAΘB,par)::Float64 = 2/π * P_primordial(par, k) * k^2 * ΘAΘB # TODO: function barrier on P_primordial?
dCl_dk_TT(l,k,S,SE,η,par)::Float64 = dCl_dk_generic(l, k, Θl0(l,k,S,η)^2,              par)
dCl_dk_TE(l,k,S,SE,η,par)::Float64 = dCl_dk_generic(l, k, Θl0(l,k,S,η)*ΘEl0(l,k,SE,η), par)
dCl_dk_EE(l,k,S,SE,η,par)::Float64 = dCl_dk_generic(l, k, ΘEl0(l,k,SE,η)^2,            par)

# TODO: loop over types, to spline only once?
function Cl(rec::Recombination, ls::Vector{Int}, type::Symbol)
    # TODO
    bg = rec.bg
    par = bg.par
    η = bg.η
    xs = range(-20, 0, length=2000) # TODO: looks like this is enough?
    logks = range(log10(1/(c*η(0))), log10(4000/(c*η(0))), length=250) # TODO: 250

    Sspl, SEspl = spline_S(rec, [S, SE], xs, logks)

    Clarr = Vector{Float64}(undef, length(ls))
    Threads.@threads for i in 1:length(ls)
        l = ls[i]

        time = @elapsed begin
        if type == :TT
            integrand = k -> dCl_dk_TT(l,k,Sspl,SEspl,η,par)
        elseif type == :TE || type == :ET
            integrand = k -> dCl_dk_TE(l,k,Sspl,SEspl,η,par)
        elseif type == :EE
            integrand = k -> dCl_dk_EE(l,k,Sspl,SEspl,η,par)
        else
            throw(error("unknown Cl type: $type"))
        end
        Clarr[i] = integrate_trapz(integrand, 1/(c*η(0)), 4000/(c*η(0)); step=2*π/(c*η(0)*10))
        #Clarr[i] = integrate_adaptive(integrand, 1/(c*η(0)), 4000/(c*η(0)))
        end
        println("Cl(l=$l) = $(Clarr[i]) ($time seconds)")
    end
    return Clarr
end

Dl(l, Cl, Tγ0) = l * (l+1) / (2*π) * Cl * (Tγ0 / 1e-6)^2 # convert to "Planck units"
function Dl(rec::Recombination, ls::Vector{Int}, type::Symbol)
    return Dl.(ls, Cl(rec, ls, type), rec.bg.par.Tγ0)
end

# TODO: spline multiple S in one pass?
function spline_S(rec::Recombination, Sfunc::Function, xs::AbstractRange, logks::AbstractRange, perturbs::Vector{PerturbationMode})
    Sdata = Float64[Sfunc(perturb, x) for x in xs, perturb in perturbs]
    Sdata[xs .== 0, :] .= 0.0 # TODO: set by interpolation instead?

    S_spline_x_logk = spline((xs, logks), Sdata) # (x, log10(k)) spline - must convert to S(x,k) spline upon calling it!
    S_spline_x_k(x, k) = S_spline_x_logk(x, log10(k))
    return S_spline_x_k
end

# Useful for splining S and SE, while computing perturbations only once
function spline_S(rec::Recombination, Sfuncs::Vector{<:Function}, xs::AbstractRange, logks::AbstractRange)
    ks = 10 .^ logks
    perturbs = [PerturbationMode(rec, k) for k in ks] # TODO: multi-threading
    return [spline_S(rec, Sfunc, xs, logks, perturbs) for Sfunc in Sfuncs]
end

function spline_S(rec::Recombination, Sfunc::Function, xs::AbstractRange, logks::AbstractRange)
    return spline_S(rec, [Sfunc], xs, logks)[1]
end

function spline_S(rec::Recombination, Sfunc)
    η = rec.bg.η
    xs = range(-20, 0, length=2000)
    logks = range(log10(1/(c*η(0))), log10(4000/(c*η(0))), length=250)
    return spline_S(rec, Sfunc, xs, logks)
end
