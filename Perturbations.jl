const polarization = true
const neutrinos = true
const lγmax = 6  # recommended by notes
const lνmax = 10 # recommended by Callin
#const lmax = max(lγmax, lνmax) # operate with only one lmax # TODO: keep lmax a function variable, then construct an object i with fields like i.δc, i.Θl(1), etc.
const lmax = 10
@assert lmax >= 4 # equations are ambiguous with lmax <= 3 (how to handle l<2, l=2 and l=lmax?

# variable index map (1-based),
# ordered so full system *extends* the tightly coupled one (useful for switching)
const i_δc = 1
const i_vc = 2
const i_δb = 3
const i_vb = 4
const i_Φ  = 5
const i_N0 = 6
const i_Nl(l::Integer) = i_N0 + l
const i_Θ0 = i_Nl(lmax) + 1
const i_Θl(l::Integer) = i_Θ0 + l
const i_max_tight = i_Θl(1)  # last variable of tight system
const i_ΘP0 = i_Θl(lmax) + 1 # variables here and below are only part of the full system
const i_ΘPl(l::Integer) = i_ΘP0 + l
const i_max = i_ΘPl(lmax)

function time_tight_coupling(co::ΛCDM, k::Real; tol::Float64=10.0)
    x1 = find_zero(x -> abs(dτ(co,x)) - tol,                (-20.0, +20.0)) # tight coupling assumes     1/τ′ << 1
    x2 = find_zero(x -> abs(dτ(co,x)) - tol * c*k/aH(co,x), (-20.0, +20.0)) # tight coupling assumes ck/̇(aτ′) << 1
    x3 = time_switch_Peebles(co) # switch no later than when recombination begins # TODO: or -8.3?
    x4 = -10.0 # switch before this to avoid kink in ΘP # TODO: find dynamic way of computing this?
    return min(x1, x2, x3, x4) # take the earliest of these times
end

function perturbations_initial_conditions(co::ΛCDM, x0::Real, k::Real)
    fν = co.Ων0 / co.Ωr0
    Ψ  = -1 / (3/2 + 2*fν/5)

    y = Vector{Float64}(undef, i_max)
    y[i_δc] = y[i_δb] = -3/2 * Ψ
    y[i_vc] = y[i_vb] = -k*c/(2*aH(co,x0)) * Ψ
    y[i_Φ]            = -(1 + 2*fν/5) * Ψ
    y[i_Θl(0)]        = -1/2 * Ψ
    y[i_Θl(1)]        = -c*k / (3*aH(co,x0)) * y[i_Θl(0)]
    y[i_Θl(2)]        = (polarization ? -8/15 : -20/45) * c*k / (aH(co,x0)*dτ(co,x0)) * y[i_Θl(1)]
    for l in 3:lmax
        y[i_Θl(l)]    = -l/(2*l+1) * c*k/(aH(co,x0)*dτ(co,x0)) * y[i_Θl(l-1)] # recursive relation
    end

    if polarization
        y[i_ΘPl(0)]       = 5/4 * y[i_Θl(2)]
        y[i_ΘPl(1)]       = -c*k/(4*aH(co,x0)*dτ(co,x0)) * y[i_Θl(2)]
        y[i_ΘPl(2)]       = 1/4 * y[i_Θl(2)]
        for l in 3:lmax
            y[i_ΘPl(l)]    = -l/(2*l+1) * c*k/(aH(co,x0)*dτ(co,x0)) * y[i_ΘPl(l-1)] # recursive relation
        end
    else
        y[i_ΘPl(0):i_ΘPl(lmax)] .= 0.0
    end

    if neutrinos
        y[i_Nl(0)] = -1/2 * Ψ
        y[i_Nl(1)] = +c*k/(6*aH(co,x0)) * Ψ
        y[i_Nl(2)] = +(c*k*a(x0)/co.H0)^2 / (30*co.Ωr0) * Ψ # expand to avoid 1/Ων0 with Ων0=0
        for l in 3:lmax
            y[i_Nl(l)] = c*k/((2*l+1)*aH(co,x0)) * y[i_Nl(l-1)]
        end
    else
        y[i_Nl(0):i_Nl(lmax)] .= 0.0
    end

    return y
end

function perturbations_mode_tight(co::ΛCDM, k::Real; x1::Real=-20.0, x2::Real=0.0, stiff=false, kwargs...)
    function dy_dx!(x, y, dy)
        # pre-compute some common combined quantities
        ck_aH = c*k / aH(co,x)
        daH_aH = daH(co,x) / aH(co,x)
        R = 4*co.Ωγ0 / (3*co.Ωb0*a(x))

        # 1) un-pack variables from vector
        δc = y[i_δc]
        vc = y[i_vc]
        δb = y[i_δb]
        vb = y[i_vb]
        Φ  = y[i_Φ]
        N0 = y[i_Nl(0)]
        Nl = @view y[i_Nl(1):i_Nl(lmax)]
        Θ0 = y[i_Θl(0)]
        Θ1 = y[i_Θl(1)]
        Θ2 = (polarization ? -8/15 : -20/45) * ck_aH * Θ1/dτ(co,x)

        # 2) compute derivatives
        Ψ   = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θ2 + co.Ων0*Nl[2])
        dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0 + 4*co.Ων0/a(x)^2*N0)
        dδc = ck_aH*vc - 3*dΦ
        dδb = ck_aH*vb - 3*dΦ
        dvc = -vc - ck_aH*Ψ
        dΘ0 = -ck_aH*Θ1 - dΦ

        # calculate (dΘ1, dΘ2) (and dvb) with fixed-point iteration, as suggested by Callin
        d_aHdτ = daH(co,x)*dτ(co,x) + aH(co,x)*d2τ(co,x) # needed in fixed-point iteration
        function dvb_dΘ1_dΘ2_fixed_point(dvb_dΘ1_dΘ2)
            _, dΘ1, dΘ2 = dvb_dΘ1_dΘ2 # ignore dvb; it's in the tuple just so it is available after the fixed-point iteration
            d_3Θ1_plus_vb = (-((1-R)*dτ(co,x)+(1+R)*d2τ(co,x))*(3*Θ1+vb) - ck_aH*Ψ + (1-daH_aH)*ck_aH*(-Θ0+2*Θ2) + ck_aH*(-dΘ0+2*dΘ2)) / # notes approximate Θ2′=0 here (see discussion under Callin (34))
                            ((1+R)*dτ(co,x) + daH_aH - 1) # Callin's 3*Θ1′ + vb′ = Hans' q
            dvb = 1/(1+R) * (-vb - ck_aH*Ψ + R*(d_3Θ1_plus_vb+ck_aH*(-Θ0+2*Θ2) - ck_aH*Ψ))
            dΘ1 = (d_3Θ1_plus_vb - dvb) / 3
            dΘ2 = (polarization ? -8/15 : -20/45) * ck_aH * (dΘ1/dτ(co,x) - Θ1*d_aHdτ / (aH(co,x)*dτ(co,x)^2)) # anal Θ2′(x) # TODO: or from q?
            return (dvb, dΘ1, dΘ2)
        end
        dvb, dΘ1, _ = fixed_point_iterate(dvb_dΘ1_dΘ2_fixed_point, (NaN, 0.0, 0.0); tol=1e-10) # don't need dΘ2 below

        # re-pack variables into vector
        dy[i_δc]    = dδc
        dy[i_vc]    = dvc
        dy[i_δb]    = dδb
        dy[i_vb]    = dvb
        dy[i_Φ]     = dΦ
        dy[i_Θl(0)] = dΘ0
        dy[i_Θl(1)] = dΘ1

        # neutrinos
        if neutrinos
            dy[i_Nl(0)] = -ck_aH*Nl[1] - dΦ
            dy[i_Nl(1)] =  ck_aH/3 * (N0 - 2*Nl[2] + Ψ)
            for l in 2:lmax-1
                dy[i_Nl(l)] = ck_aH/(2*l+1) * (l*Nl[l-1] - (l+1)*Nl[l+1])
            end
            dy[i_Nl(lmax)] = ck_aH*Nl[lmax-1] - (lmax+1)/(aH(co,x)*η(co,x)) * Nl[lmax]
        else
            dy[i_Nl(0):i_Nl(lmax)] .= 0.0
        end

        return nothing
    end

    splines = Vector{Spline1D}(undef, i_max)
    alg_hints = [stiff ? :auto : :nonstiff]
    y1 = perturbations_initial_conditions(co, x1, k)[1:i_max_tight] # cut away variables not in the tight system
    splines[1:i_max_tight] .= _spline_integral(dy_dx!, x1, x2, y1; alg_hints=alg_hints, name="perturbations tight (k=$(k*Mpc)/Mpc)", kwargs...)

    # extend tight splines to full system
    x = splinex(splines[1])
    splines[i_Θl(2)] = Spline1D(x, @. (polarization ? -8/15 : -20/45) * c*k / (aH(co,x)*dτ(co,x)) * splines[i_Θl(1)](x); bc="error")
    for l in 3:lmax
        splines[i_Θl(l)] = Spline1D(x, @. -l/(2*l+1) * c*k/(aH(co,x)*dτ(co,x)) * splines[i_Θl(l-1)](x); bc="error")
    end

    if polarization
        splines[i_ΘPl(0)] = Spline1D(x, @. 5/4 * splines[i_Θl(2)](x); bc="error")
        splines[i_ΘPl(1)] = Spline1D(x, @. -c*k/(4*aH(co,x)*dτ(co,x)) * splines[i_Θl(2)](x); bc="error")
        splines[i_ΘPl(2)] = Spline1D(x, @. 1/4 * splines[i_Θl(2)](x); bc="error")
        for l in 3:lmax
            splines[i_ΘPl(l)] = Spline1D(x, @. -l/(2*l+1) * c*k/(aH(co,x)*dτ(co,x)) * splines[i_ΘPl(l-1)](x); bc="error")
        end
    else
        for l in 0:lmax
            splines[i_ΘPl(l)] = Spline1D(x, 0 .* x)
        end
    end

    return splines
end

function perturbations_mode_full(co::ΛCDM, k::Real; y1=nothing, x1::Real=-20.0, x2::Real=0.0, stiff=true, kwargs...)
    function dy_dx!(x::Float64, y::Vector{Float64}, dy::Vector{Float64})
        # pre-compute some common combined quantities
        ck_aH = c*k / aH(co,x)
        R = 4*co.Ωγ0 / (3*co.Ωb0*a(x))
        τ′ = dτ(co,x) # ′ = \prime != '

        # 1) un-pack variables from vector
        δc = y[i_δc]
        vc = y[i_vc]
        δb = y[i_δb]
        vb = y[i_vb]
        Φ  = y[i_Φ]
        N0 = y[i_Nl(0)]
        Nl = @view y[i_Nl(1):i_Nl(lmax)]
        Θ0 = y[i_Θl(0)]
        Θl = @view y[i_Θl(1):i_Θl(lmax)] # for l ≥ 1 (since Julia is 1-indexed) # TODO: use OffsetArrays?
        ΘP0 = y[i_ΘPl(0)]
        ΘPl = @view y[i_ΘPl(1):i_ΘPl(lmax)]

        # 2) compute derivatives
        Ψ   = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θl[2]+co.Ων0*Nl[2])
        Π   = Θl[2] + ΘP0 + ΘPl[2]
        dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0 + 4*co.Ων0/a(x)^2*N0)
        dδc = ck_aH*vc - 3*dΦ
        dδb = ck_aH*vb - 3*dΦ
        dvc = -vc - ck_aH*Ψ
        dvb = -vb - ck_aH*Ψ + τ′*R*(3*Θl[1]+vb)
        dΘl0   = -ck_aH*Θl[1] - dΦ
        dΘl1   =  ck_aH/3 * (Θ0-2*Θl[2]+Ψ) + τ′*(Θl[1]+vb/3)
        dΘlmax = ck_aH*Θl[lmax-1] - (lmax+1)/(aH(co,x)*η(co,x))*Θl[lmax] + τ′*Θl[lmax] # 2nd term: their η is my c*η
        dΘPl0   = -ck_aH*ΘPl[1] + τ′*(ΘP0-Π/2)
        dΘPlmax =  ck_aH*ΘPl[lmax-1] - (lmax+1)/(aH(co,x)*η(co,x))*ΘPl[lmax] + τ′*ΘPl[lmax]

        # re-pack variables into vector
        dy[i_δc]    = dδc
        dy[i_vc]    = dvc
        dy[i_δb]    = dδb
        dy[i_vb]    = dvb
        dy[i_Φ]     = dΦ
        dy[i_Θl(0)] = dΘl0
        dy[i_Θl(1)] = dΘl1
        for l in 2:lmax-1
            dy[i_Θl(l)] = ck_aH/(2*l+1) * (l*Θl[l-1] - (l+1)*Θl[l+1]) + τ′*(Θl[l]-Π/10*δ(l,2))
        end
        dy[i_Θl(lmax)] = dΘlmax

        if polarization
            dy[i_ΘPl(0)] = dΘPl0
            dy[i_ΘPl(1)] = ck_aH/(2*1+1) * (1*ΘP0 - (1+1)*ΘPl[1+1]) + τ′*(ΘPl[1]-Π/10*δ(1,2))
            for l in 2:lmax-1
                dy[i_ΘPl(l)] = ck_aH/(2*l+1) * (l*ΘPl[l-1] - (l+1)*ΘPl[l+1]) + τ′*(ΘPl[l]-Π/10*δ(l,2))
            end
            dy[i_ΘPl(lmax)] = dΘPlmax
        else
            dy[i_ΘPl(0):i_ΘPl(lmax)] .= 0.0
        end

        # neutrinos
        if neutrinos
            dy[i_Nl(0)] = -ck_aH*Nl[1] - dΦ
            dy[i_Nl(1)] =  ck_aH/3 * (N0 - 2*Nl[2] + Ψ)
            for l in 2:lmax-1
                dy[i_Nl(l)] = ck_aH/(2*l+1) * (l*Nl[l-1] - (l+1)*Nl[l+1])
            end
            dy[i_Nl(lmax)] = ck_aH*Nl[lmax-1] - (lmax+1)/(aH(co,x)*η(co,x)) * Nl[lmax]
        else
            dy[i_Nl(0):i_Nl(lmax)] .= 0.0
        end

        return nothing # dy is in-place
    end

    if isnothing(y1)
        y1 = perturbations_initial_conditions(co, x1, k)
    end
    alg_hints = [stiff ? :auto : :nonstiff]
    return _spline_integral(dy_dx!, x1, x2, y1; alg_hints=alg_hints, name="perturbations full (k=$(k*Mpc)/Mpc)", kwargs...)
end

function perturbations_mode(co::ΛCDM, k::Real; tight::Bool=false)
    if tight
        # 1) use tight coupling approximation at early times to avoid stiff equations,
        # 2) then integrate the full equations from when the approximation breaks down,
        # 3) merge 1+2
        x12 = time_tight_coupling(co, k)
        spl1s = perturbations_mode_tight(co, k; x2=x12, stiff=false, abstol=1e-9, reltol=1e-9)
        y12 = [spl1(x12) for spl1 in spl1s] # give final tight values as ICs for full system
        spl2s = perturbations_mode_full(co, k; x1=x12, y1=y12, stiff=false, abstol=1e-9, reltol=1e-9)
        return [splinejoin(spl1s[i], spl2s[i]) for i in 1:i_max] # join splines
    else
        # only integrate the full (stiff) equations using an appropriate solver
        # TODO: use lower tolerance at later times (when system is not stiff) to speed up?
        return perturbations_mode_full(co, k; stiff=true, abstol=1e-9, reltol=1e-9)
    end
end

function perturbations_splines(co::ΛCDM; tight::Bool=false)
    if length(co.perturbation_splines) == 0
        t1 = now()

        kmin, kmax = 0.00005 / Mpc, 0.3 / Mpc
        ks = kmin .+ (kmax-kmin) * range(0, 1; length=200) .^ 2 # quadratic spacing
        
        # (x,k) spline requires 1D arrays for x and k (and does not accept a fully irregular 2D grid)
        # use the x-values from the mode that has the most (not necessarily k=kmax)
        spliness = Vector{Vector{Spline1D}}(undef, length(ks))
        @threads for i_k in 1:length(ks) # call julia with --threads=8 to get a decent speed-up
            k = ks[i_k]
            spliness[i_k] = perturbations_mode(co, k; tight=tight)
        end
        xs = splinex(spliness[argmax(length(splinex(splines[1])) for splines in spliness)][1])
        #xs = xs[1:4:end] # TODO: reduce memory usage here or in perturbation_mode?

        perturbs = Array{Float64, 3}(undef, i_max, length(xs), length(ks)) # indexed as [i_quantity, i_x, i_k]
        for i_k in 1:length(ks)
            for i_qty in 1:i_max
                perturbs[i_qty, :, i_k] .= spliness[i_k][i_qty](xs)
            end
        end

        for i_qty in 1:i_max
            qty = @view perturbs[i_qty, :, :]
            push!(co.perturbation_splines, Spline2D(xs, ks, qty)) # spline (x, k)
        end

        t2 = now()
        println("Splined perturbations(x,k) on ($(length(xs)),$(length(ks))) grid in $(t2-t1) using $(sizeof(perturbs)/1e6) MB and $(nthreads()) parallel threads")
    end
    return co.perturbation_splines
end

Φ(co::ΛCDM, x::Real, k::Real) = perturbations_splines(co)[i_Φ](x, k)
Ψ(co::ΛCDM, x::Real, k::Real) = -Φ(co,x,k) - 12*co.H0^2/(c*k*a(x))^2 * (co.Ωγ0*Θl(co,x,k,2) + co.Ων0*Nl(co,x,k,2))
δc(co::ΛCDM, x::Real, k::Real) = perturbations_splines(co)[i_δc](x, k)
δb(co::ΛCDM, x::Real, k::Real) = perturbations_splines(co)[i_δb](x, k)
δγ(co::ΛCDM, x::Real, k::Real) = 0 # TODO
vc(co::ΛCDM, x::Real, k::Real) = perturbations_splines(co)[i_vc](x, k)
vb(co::ΛCDM, x::Real, k::Real) = perturbations_splines(co)[i_vb](x, k)
Θl( co::ΛCDM, x::Real, k::Real, l::Integer) = perturbations_splines(co)[i_Θl(l)](x, k)
ΘPl(co::ΛCDM, x::Real, k::Real, l::Integer) = perturbations_splines(co)[i_ΘPl(l)](x, k)
Nl( co::ΛCDM, x::Real, k::Real, l::Integer) = perturbations_splines(co)[i_Nl(l)](x, k)
