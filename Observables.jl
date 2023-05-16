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

#=
function Cl(l::Integer, η, par, Θ0A, Θ0B; integratek=integrate_adaptive)
    dCl_dk(k)::Float64 = 2/π * P_primordial(par, k) * k^2 * Θ0A(k) * Θ0B(k) # TODO: function barrier on P_primordial?
    Cl = integratek(dCl_dk, 1/(c*η(0)), 4000/(c*η(0)))
    return Cl
end

function Cl(rec::Recombination, ls::Vector{Int}, Θ0A, Θ0B)
    println("Entered Cl")
    Clarr = Vector{Float64}(undef, length(ls))
    #Threads.@threads for i in 1:length(ls)
    for i in 1:length(ls)
        l = ls[i]
        println("C(l=$l) = ")
        Θ0Af(k::Float64)::Float64 = Θ0A(l,k)
        Θ0Bf(k::Float64)::Float64 = Θ0B(l,k)
        Clarr[i] = Cl(l, rec.bg.η, rec.bg.par, Θ0Af, Θ0Bf)
        #time = @elapsed(Clarr[i] = Cl(l, rec.bg.η, rec.bg.par, k -> Θ0A(l,k), k -> Θ0A(l,k)))
        #println("$(Clarr[i]) ($time seconds)")
    end
    return Clarr
end

function ClTT(rec::Recombination, η::ODESolution, ls::Vector{Int})
    xs = range(-20, 0, length=2000) # TODO: looks like this is enough?
    logks = range(log10(1/(c*η(0))), log10(4000/(c*η(0))), length=250) # TODO: 250
    Sspl = spline_S(rec, S, xs, logks) # (x,k)-callable
    println("Splined S")
    Θ0A(l,k) = Θl0(l, k, Sspl, η)^2
    Θ0B(l,k) = 1.0
    return Cl(rec, ls, Θ0A, Θ0B)
end

ClTT(rec::Recombination, ls::Vector{Int}) = ClTT(rec, rec.bg.η, ls)
=#

#=
function ClTE(l::Integer, SE, η, par)
    η0 = η(0.0)
    integratex = integrate_adaptive
    dΘ0_dx(x, k) :: Float64 = SE(x, k) * jl(l, c*k*(η0-η(x)))
    Θ0(k) = √((l+2)*(l+1)*l*(l-1)) * integratex(x -> dΘ0_dx(x, k), -20.0, 0.0)
    return Cl(l, S, η, par, k -> Θ0(k)^2, k -> 1.0)
end
=#

# TODO: would be quicker without two overloaded functions
#=
function Cl(rec::Recombination, ls::Vector{Int}, Θ0A, Θ0B; integrate=integratek, verbose=false)
    bg = rec.bg
    par = bg.par
    η = bg.η

    Clarr = Vector{Float64}(undef, length(ls))
    #Threads.@threads for i in 1:length(ls)
    for i in 1:length(ls)
        l = ls[i]

        time = @elapsed begin
        dCl_dk(k)::Float64 = 2/π * P_primordial(par, k) * k^2 * Θ0A(l,k) * Θ0B(l,k) # TODO: function barrier on P_primordial?
        Clarr[i] = integrate(dCl_dk, 1/(c*η(0)), 4000/(c*η(0)))
        end
        if verbose
            println("Cl(l=$l) = $(Clarr[i]) ($time seconds)")
        end
    end
    return Clarr
end
=#

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

#=
function Cl_trapz(l, Sspl, η, par, Ny_per_osc=20, Nk_per_osc=10)
    η0 = η(0.0)

    # TODO: allocate x-array ONCE

    Ny = 1000
    dΘ0_dk_grid = Vector{Float64}(undef, Ny)
    function Θ0(k)
        # TODO: set ymin, ymax instead?
        #xmin = -20.0
        #xmax = 0.0
        #y(x) = c*k*(η0-η(x))
        #ymin = y(xmax)
        #ymax = y(xmin)
        #Δy = 2*π/Ny_per_osc
        #Ny = Int(round((ymax-ymin) / Δy))
        #Ny = max(Ny, 1000) # TODO: do the maximum with what is needed to resolve source function
        #Ny = 1000
        xs = range(xmin, xmax, length=Ny)
        ys = y.(xs)
        dΘ0_dk_grid .= Sspl.(xs, k) .* jl.(l, c*k * (η0 .- η.(xs))) # (x,k) function
        return trapz(xs, dΘ0_dk_grid) # integrate over x
    end

    # TODO: spline jl. all allocations are from it!

    # TODO: @code_warntype gives η(xs)::Any ???!?!? use η.(xs) instead?
    # TODO: jl allocates a lot. reduce somehow?
    #dΘ0_dk_grid = evalgrid(Sspl, xs, ks) .* jl.(l, c * ks' .* (η(0.0) .- η.(xs))) # (x,k) function

    kmin =    1 / (c*η(0))
    kmax = 4000 / (c*η(0))
    Δk = 2*π / (c*η0*Nk_per_osc)
    ks = range(kmin, kmax, step=Δk)

    println("$(length(ks)) k-values")

    dCl_dk_grid = 2/π * P_primordial.(par, ks) .* (ks .* Θ0.(ks)) .^ 2
    return trapz(ks, dCl_dk_grid) # integrate over k
end
=#

# TODO: for some reason it fails for l=8 and rtol=1e-1
#=
function Cl_cubature(l, Sspl, η, par; rtol=1e-3, kwargs...)
    η0 = η(0.0)
    S(x,k) = Sspl(x, log10(k)) # TODO: log or not?
    function integrand(vec)
        k, x1, x2 = vec[1], vec[2], vec[3]
        return S(x1, k) * jl(l, c*k*(η0 - η(x1))) * # dΘ0_dx (#1)
               S(x2, k) * jl(l, c*k*(η0 - η(x2))) * # dΘ0_dx (#2)
               2/π * k^2 * P_primordial(par, k)   # dCl_dk
    end
    kmin = 1 / (c*η0)
    kmax = 4000 / (c*η0)
    xmin = -20.0
    xmax = 0.0
    return hcubature(integrand, (kmin, xmin, xmin), (kmax, xmax, xmax); rtol=rtol, kwargs...)[1]
end
=#

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
