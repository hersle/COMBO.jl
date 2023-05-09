P_primordial(par::Parameters, k::Real) = 2*π^2/k^3 * par.As * (k/par.k_pivot)^(par.ns-1)
Δ(perturb::Perturbations, x::Real, k::Real) = c^2*k^2*Φ(perturb,x) / (3/2*perturb.rec.bg.par.Ωm0/a(x)*perturb.rec.bg.par.H0^2)
P(perturb::Perturbations, x::Real, k::Real) = abs(Δ(perturb,x,k))^2 * P_primordial(perturb.rec.bg.par, k)

function Cls(rec::Recombination, ls::Vector)
    bg = rec.bg

    kmin =    1 / (c*η(bg,0))
    kmax = 3000 / (c*η(bg,0))
    Δk = 2*π / (c*η(bg,0)*20)
    ks = range(kmin, kmax, step=Δk)
    iks = 1:length(ks)

    perturbs = [Perturbations(rec, k) for k in ks]

    # TODO: this is extremely slow now. put splines behind function barriers
    x = range(-20, 0, length=1000) # TODO: take from perturbation solution?
    Θl0(l, ik) = trapz(x, S.(perturbs[ik],x) .* jl.(l, c*ks[ik]*(η(bg,0.0) .- η.(bg,x)))) # TODO: or trapz?

    Clarr = []
    for l in ls
        println("l = $l")
        push!(Clarr, Cl(l))
    end

    return Clarr
    #return Cl.(ls)

    #=
    println("l = $l")
    kmin =    1 / (c*η(bg,0))
    kmax = 3000 / (c*η(bg,0))
    Δk = 2*π / (c*η(bg,0)*10)
    k = range(kmin, kmax, step=Δk)
    x = range(-20, 0, length=1000)
    Θl0(l, k) = trapz(x, S(perturb,x) .* jl.(l, c*k*(η(bg,0.0) .- η.(bg,x)))) # TODO: or trapz?
    return 2/π * trapz(k, @. k^2 * P_primordial(co,k) * Θl0(l,k)^2)
    =#
end
