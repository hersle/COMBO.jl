# TODO: power spectrum struct?
export P_primordial, P, Cl, Dl

using Interpolations # splines
using Bessels: sphericalbesselj as jl # spherical Bessel function # TODO: correct function?
using Base.Threads # parallelization
using Trapz # quadrature

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

# TODO: use dynamic x grid depending on k and l?
dΘTl0_dx(l, xs, ks, STs, η) =                              STs .* jl.(l, c * ks' .* (η(0) .- η.(xs))) # TODO: pass array of η?
dΘEl0_dx(l, xs, ks, SEs, η) = √((l+2)*(l+1)*(l+0)*(l-1)) * SEs .* jl.(l, c * ks' .* (η(0) .- η.(xs)))

ΘTl0(l, xs, ks, STs, η) = trapz(xs, dΘTl0_dx(l,xs,ks,STs,η), Val(1)) # integrate over x (axis #1), not k (axis #2)
ΘEl0(l, xs, ks, SEs, η) = trapz(xs, dΘEl0_dx(l,xs,ks,SEs,η), Val(1)) # integrate over x (axis #1), not k (axis #2)

dCl_dk_generic(l,ks,ΘAΘBs,par) = 2/π * P_primordial.(par, ks) .* ks .^ 2 .* ΘAΘBs # TODO: function barrier on P_primordial?
dCl_dk_TT(l,xs,ks,STs,SEs,η,par) = dCl_dk_generic(l, ks, ΘTl0(l,xs,ks,STs,η) .^ 2,                   par)
dCl_dk_TE(l,xs,ks,STs,SEs,η,par) = dCl_dk_generic(l, ks, ΘTl0(l,xs,ks,STs,η) .* ΘEl0(l,xs,ks,SEs,η), par)
dCl_dk_EE(l,xs,ks,STs,SEs,η,par) = dCl_dk_generic(l, ks, ΘEl0(l,xs,ks,SEs,η) .^ 2,                   par)

trapz_extra(x0, y0, x, y) = trapz(x, y) + (x[1]-x0) * (y0 + y[1]) / 2 # trapezoid integral of (x, y) extended with another leftmost point (x0, y0)
Cl_TT(l,xs,ks,STs,SEs,η,par) = trapz_extra(0.0, 0.0, ks, dCl_dk_TT(l,xs,ks,STs,SEs,η,par)) # integrate over k (manually add k=0)
Cl_TE(l,xs,ks,STs,SEs,η,par) = trapz_extra(0.0, 0.0, ks, dCl_dk_TE(l,xs,ks,STs,SEs,η,par)) # integrate over k (manually add k=0)
Cl_EE(l,xs,ks,STs,SEs,η,par) = trapz_extra(0.0, 0.0, ks, dCl_dk_EE(l,xs,ks,STs,SEs,η,par)) # integrate over k (manually add k=0)

# TODO: loop over types, to spline only once?
function Cl(rec::Recombination, ls::Vector{Int}, type::Symbol; spline_S_before_gridding=true)
    # TODO
    bg = rec.bg
    par = bg.par
    η = bg.η
    xs = range(-10, 0, step=0.01) # TODO: looks like this is enough?
    ks = range(1/(c*η(0)), 4000/(c*η(0)), step=2*π/(c*η(0)*10))

    STs, SEs = grid_S(rec, [S, SE], xs, ks; spline_first=spline_S_before_gridding)

    Clarr = Vector{Float64}(undef, length(ls))
    Threads.@threads for i in 1:length(ls)
        l = ls[i]

        time = @elapsed begin
        if type == :TT
            Clarr[i] = Cl_TT(l, xs, ks, STs, SEs, η, par)
        elseif type == :TE
            Clarr[i] = Cl_TE(l, xs, ks, STs, SEs, η, par)
        elseif type == :EE
            Clarr[i] = Cl_EE(l, xs, ks, STs, SEs, η, par)
        else
            throw(error("unknown Cl type: $type"))
        end
        end
        println("Cl(l=$l) = $(Clarr[i]) ($time seconds)")
    end
    return Clarr
end

Dl(l, Cl, Tγ0) = l * (l+1) / (2*π) * Cl * Tγ0^2 # convert to "Planck units"
function Dl(rec::Recombination, ls::Vector{Int}, type::Symbol)
    return Dl.(ls, Cl(rec, ls, type), rec.bg.par.Tγ0)
end

# Useful for splining S and SE, while computing perturbations only once
function spline_S(rec::Recombination, Sfuncs::Vector{<:Function}, xs::AbstractRange, logks::AbstractRange)
    Sgrids = grid_S(rec, Sfuncs, xs, 10 .^ logks)
    Sspls = [spline((xs, logks), Sgrid) for Sgrid in Sgrids] # (x, log(k)) callables
    return [(x, k) -> Sspl(x, log10(k)) for Sspl in Sspls] # (x, k) callables (this is what the "user" wants!)
end

function spline_S(rec::Recombination, Sfunc::Function, xs::AbstractRange, logks::AbstractRange)
    return spline_S(rec, [Sfunc], xs, logks)[1]
end

function spline_S(rec::Recombination, Sfunc)
    η = rec.bg.η
    xs = range(-20, 0, step=0.01)
    logks = range(log10(1/(c*η(0))), log10(4000/(c*η(0))), length=300)
    return spline_S(rec, Sfunc, xs, logks)
end

function grid_S(rec::Recombination, Sfuncs::Vector{<:Function}, xs, ks; spline_first=false) :: Vector{Matrix{Float64}}
    if spline_first
        Sspls = spline_S(rec, Sfuncs)
        return [Sspl.(xs, ks') for Sspl in Sspls]
    else
        perturbs = Vector{PerturbationMode}(undef, length(ks))
        Threads.@threads for i in 1:length(ks)
            perturbs[i] = PerturbationMode(rec, ks[i])
        end
        # TODO: force S(x=0) = 0
        return [[Sfunc(perturb, x) for x in xs, perturb in perturbs] for Sfunc in Sfuncs]
    end
end

function grid_S(rec::Recombination, Sfunc::Function, xs, ks; spline_first=false)
    return grid_S(rec, [Sfunc], xs, ks; spline_first=spline_first)[1]
end
