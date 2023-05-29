using Interpolations # splines
using Bessels: sphericalbesselj as jl # spherical first-kind Bessel function
using Base.Threads # parallelization
using Trapz # quadrature

export MatterPowerSpectrum
export P_primordial, P

export CMBPowerSpectrum
export Cl, Cl_TT, Cl_EE, Cl_TE, Dl, dCl_dk_TT, dCl_dk_TE, dCl_dk_EE
export ΘTl0, ΘEl0, dΘTl0_dx, dΘEl0_dx

struct MatterPowerSpectrum
    rec::Recombination
    k::Vector{Float64}
    Pk::Vector{Float64}

    function MatterPowerSpectrum(rec::Recombination, ks::Vector{Float64}=10 .^ range(-4,+2,length=200) / Mpc)
        perturbs = PerturbationModes(rec, ks)
        Pks = [P(perturb, 0) for perturb in perturbs]
        new(rec, ks, Pks)
    end
end

P_primordial(par::Parameters, k) = 2*π^2/k^3 * par.As * (k/par.k_pivot)^(par.ns-1)
Δ(perturb::PerturbationMode, x) = c^2*perturb.k^2*Φ(perturb,x) / (3/2*perturb.rec.bg.par.Ωm0/a(x)*perturb.rec.bg.par.H0^2)
P(perturb::PerturbationMode, x) = abs(Δ(perturb,x))^2 * P_primordial(perturb.rec.bg.par, perturb.k)

struct CMBPowerSpectrum
    rec::Recombination
    l::Vector{Int}
    Cl::Vector{Float64}
    Dl::Vector{Float64}

    function CMBPowerSpectrum(rec::Recombination, ls::Vector{Int}, type::Symbol; xs=range(-10, 0, step=0.01), kcη0s=range(1, 4000, step=2*π/10), spline_S_before_gridding=true, verbose=true)
        @assert type in (:TT, :EE, :TE) "unknown Cl type: $type"
        bg = rec.bg
        par = bg.par
        ks = kcη0s / (c*η(bg,0))
        χs = χ.(bg, xs) # = c * (η(0) - η(x))
        STs, SEs = grid_S(rec, [ST, SE], xs, ks; spline_first=spline_S_before_gridding)

        Cls = Vector{Float64}(undef, length(ls))
        Threads.@threads for i in 1:length(ls)
            l = ls[i]
            time = @elapsed begin

            if type == :TT
                Cls[i] = Cl_TT(l, xs, ks, STs, SEs, χs, par)
            elseif type == :TE
                Cls[i] = Cl_TE(l, xs, ks, STs, SEs, χs, par)
            elseif type == :EE
                Cls[i] = Cl_EE(l, xs, ks, STs, SEs, χs, par)
            end
            end
            if verbose
                println("Cl(l=$l) = $(Cls[i]) ($time seconds)")
            end
        end

        Dls = Dl.(ls, Cls, rec.bg.par.Tγ0)
        new(rec, ls, Cls, Dls)
    end

    function CMBPowerSpectrum(rec::Recombination, type::Symbol; n1=10, n2=20, n3=150, kwargs...)
        ls = unique(Int.(round.(10 .^ vcat(range(0.0, 1.0, length=n1),
                                           range(1.0, 2.0, length=n2),
                                           range(2.0, 3.4, length=n3)))))
        return CMBPowerSpectrum(rec, ls, type; kwargs...)
    end
end

# TODO: 1) use dynamic x grid depending on k and l?
# TODO: 2) use Levin integration (see https://discourse.julialang.org/t/rfc-ann-oscillatoryintegralsode-jl-levin-method-ordinarydiffeq/55601)?
dΘTl0_dx(l, xs, ks, STs, χs) =                              STs .* jl.(l, ks' .* χs)
dΘEl0_dx(l, xs, ks, SEs, χs) = √((l+2)*(l+1)*(l+0)*(l-1)) * SEs .* jl.(l, ks' .* χs)

ΘTl0(l, xs, ks, STs, χs) = trapz(xs, dΘTl0_dx(l,xs,ks,STs,χs), Val(1)) # integrate over x (axis #1), not k (axis #2)
ΘEl0(l, xs, ks, SEs, χs) = trapz(xs, dΘEl0_dx(l,xs,ks,SEs,χs), Val(1))

dCl_dk_generic(l,ks,ΘAΘBs,par) = 2/π * P_primordial.(par, ks) .* ks .^ 2 .* ΘAΘBs
dCl_dk_TT(l,xs,ks,STs,SEs,χs,par) = dCl_dk_generic(l, ks, ΘTl0(l,xs,ks,STs,χs) .^ 2,                   par)
dCl_dk_TE(l,xs,ks,STs,SEs,χs,par) = dCl_dk_generic(l, ks, ΘTl0(l,xs,ks,STs,χs) .* ΘEl0(l,xs,ks,SEs,χs), par)
dCl_dk_EE(l,xs,ks,STs,SEs,χs,par) = dCl_dk_generic(l, ks, ΘEl0(l,xs,ks,SEs,χs) .^ 2,                   par)

trapz_extra(x0, y0, x, y) = trapz(x, y) + (x[1]-x0) * (y0 + y[1]) / 2 # trapezoid integral of (x, y) extended with another leftmost point (x0, y0)
Cl_TT(l,xs,ks,STs,SEs,χs,par) = trapz_extra(0.0, 0.0, ks, dCl_dk_TT(l,xs,ks,STs,SEs,χs,par)) # integrate over k (manually add k=0)
Cl_TE(l,xs,ks,STs,SEs,χs,par) = trapz_extra(0.0, 0.0, ks, dCl_dk_TE(l,xs,ks,STs,SEs,χs,par)) # integrate over k (manually add k=0)
Cl_EE(l,xs,ks,STs,SEs,χs,par) = trapz_extra(0.0, 0.0, ks, dCl_dk_EE(l,xs,ks,STs,SEs,χs,par)) # integrate over k (manually add k=0)

Dl(l, Cl, Tγ0) = l * (l+1) / (2*π) * Cl * Tγ0^2 # convert to "Planck units"

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
        perturbs = PerturbationModes(rec, ks)
        return [[Sfunc(perturb, x) for x in xs, perturb in perturbs] for Sfunc in Sfuncs]
    end
end

function grid_S(rec::Recombination, Sfunc::Function, xs, ks; spline_first=false)
    return grid_S(rec, [Sfunc], xs, ks; spline_first=spline_first)[1]
end
