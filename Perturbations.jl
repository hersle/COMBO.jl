# TODO: assert lmax >= 2 in main perturbations function (?) (always use this?)

# index map (1-based)
# TODO: order so one can set some y_untight from y_tight
const i_δc = 1
const i_vc = 2
const i_δb = 3
const i_vb = 4
const i_Φ  = 5
const i_Θ0 = 6
const i_Θl(l::Integer) = i_Θ0 + l
const i_max(lmax::Integer) = i_Θl(lmax)

# TODO: include neutrinos
# TODO: include polarization
# TODO: check units
# TODO: set all Θl, remove _tight suffix, then just remove them during tight coupling
function initial_conditions_tight(co::ΛCDM, x0::Real, k::Real)
    aH0 = aH(co, x0)
    dτ0 = dτ(co, x0)

    fν = 0 # TODO: neutrinos: co.Ων0 / co.Ωr0
    Ψ  = -1 / (3/2 + 2*fν/5)
    Φ  = -(1 + 2*fν/5) * Ψ # TODO: "acts as normalization" ?
    δb = δc = -3/2 * Ψ
    vb = vc = -k*c/(2*aH0) * Ψ # TODO: units?
    Θ0 = -1/2 * Ψ
    Θ1 = -c*k / (3*aH0) * Θ0

    # integration always begins with tight coupling, during which Θ(l>1) follows directly from Θ(l<=1)
    lmax = 1
    y = Vector{Float64}(undef, i_max(lmax))
    y[i_δc]    = δc
    y[i_vc]    = vc
    y[i_δb]    = δb
    y[i_vb]    = vb
    y[i_Φ]     = Φ
    y[i_Θl(0)] = Θ0
    y[i_Θl(1)] = Θ1
    return y
end

function time_tight_coupling(co::ΛCDM, k::Real)
    x1 = find_zero(x -> abs(dτ(co,x)) - 10,                (-20, +20))
    x2 = find_zero(x -> abs(dτ(co,x)) - 10 * c*k/aH(co,x), (-20, +20))
    x3 = co.x_tight_latest # TODO: -8.3? use something like time_recombination(co) or time_switch_Peebles(co) instead?
    return min(x1, x2, x3)
end

# tight coupling is equivalent to lmax=2, then post-compute l>lmax
# TODO: store splines for each requested k-value?
function splined_perturbations_tight(co::ΛCDM, k::Real; x1::Real=-20.0, x2::Real=time_tight_coupling(co, k), lmax::Integer=30)
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
        Θ0 = y[i_Θl(0)]
        Θ1 = y[i_Θl(1)]
        Θ2 = -20/45 * ck_aH * Θ1/dτ(co,x) # TODO: polarization

        # 2) compute derivatives
        Ψ   = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θ2) # TODO: neutrinos
        dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0) # TODO: neutrinos
        dδc = ck_aH*vc - 3*dΦ
        dδb = ck_aH*vb - 3*dΦ
        dvc = -vc - ck_aH*Ψ
        dΘ0 = -ck_aH*Θ1 - dΦ

        # calculate (dΘ1, dΘ2) (and dvb) with fixed-point iteration
        d_aHdτ = daH(co,x)*dτ(co,x) + aH(co,x)*d2τ(co,x) # needed in fixed-point iteration
        function dvb_dΘ1_dΘ2_fixed_point(dvb_dΘ1_dΘ2)
            _, dΘ1, dΘ2 = dvb_dΘ1_dΘ2 # ignore dvb; it's in the tuple just so it is available after the fixed-point iteration
            d_3Θ1_plus_vb = (-((1-R)*dτ(co,x)+(1+R)*d2τ(co,x))*(3*Θ1+vb) - ck_aH*Ψ + (1-daH_aH)*ck_aH*(-Θ0+2*Θ2) + ck_aH*(-dΘ0+2*dΘ2)) / # approximate Θ2' = 0 (see discussion under Callin (34))
                            ((1+R)*dτ(co,x) + daH_aH - 1) # Callin's 3*Θ1′ + vb′ = Hans' q
            dvb = 1/(1+R) * (-vb - ck_aH*Ψ + R*(d_3Θ1_plus_vb+ck_aH*(-Θ0+2*Θ2) - ck_aH*Ψ))
            dΘ1 = (d_3Θ1_plus_vb - dvb) / 3
            dΘ2 = -20/45*ck_aH * (dΘ1/dτ(co,x) - Θ1*d_aHdτ / (aH(co,x)*dτ(co,x)^2)) # anal Θ2′(x) # TODO: or from q?
            return (dvb, dΘ1, dΘ2)
        end
        dvb, dΘ1, dΘ2 = fixed_point_iterate(dvb_dΘ1_dΘ2_fixed_point, (NaN, 0.0, 0.0); tol=1e-30)

        # re-pack variables into vector
        #lmax = 1 # tight coupling: only need to integrate Θ(l<=2) (see above comment)
        dy[i_δc]    = dδc
        dy[i_vc]    = dvc
        dy[i_δb]    = dδb
        dy[i_vb]    = dvb
        dy[i_Φ]     = dΦ
        dy[i_Θl(0)] = dΘ0
        dy[i_Θl(1)] = dΘ1
        return nothing
    end

    # TODO: spline perturb(x, k)
    y1 = initial_conditions_tight(co, x1, k)
    splines = _spline_integral(dy_dx!, x1, x2, y1; abstol=1e-9, reltol=1e-9, name="perturbations tight (k=$(k*Mpc)/Mpc)")
    x = splinex(splines[1])

    # extend Θl splines up to lmax
    splines_ext = Vector{Spline1D}(undef, i_max(lmax))
    for i in 1:i_max(1)
        splines_ext[i] = splines[i]
    end
    splines_ext[i_Θl(2)] = Spline1D(x, @. -20*c*k / (45*aH(co,x)*dτ(co,x)) * splines_ext[i_Θl(1)](x)) # TODO: polarization
    for l in 3:lmax
        splines_ext[i_Θl(l)] = Spline1D(x, @. -l/(2*l+1) * c*k/(aH(co,x)*dτ(co,x)) * splines_ext[i_Θl(l-1)](x)) # recursive relation
    end
    return splines_ext
end

function initial_conditions_untight(co::ΛCDM, x0::Real, k::Real, lmax::Integer)
    # TODO: order indices so this can simply "return perturbations_tight(co, x0, k), then set the additional l>1
    return [spline(x0) for spline in splined_perturbations_tight(co, k)]
end

function splined_perturbations_untight(co::ΛCDM, k::Real; x2::Real=0.0, lmax::Integer=30) :: Vector{Spline1D}# lmax ≈ 30 (https://arxiv.org/pdf/1104.2933.pdf)
    x1 = time_tight_coupling(co, k)

    function dy_dx!(x::Float64, y::Vector{Float64}, dy::Vector{Float64})
        #print("x = $x\r") # print progress
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
        Θ0 = y[i_Θl(0)]
        Θl = @view y[i_Θl(1):i_Θl(lmax)] # for l ≥ 1 (since Julia is 1-indexed) # TODO: use OffsetArrays?

        # 2) compute derivatives
        Ψ   = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θl[2]) # TODO: Ωγ0 (Hans) or Ωr (Callin)? TODO: neutrinos
        Π   = Θl[2] + 0 # TODO: include or not? # TODO: polarization
        dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0) # TODO: neutrinos
        dδc = ck_aH*vc - 3*dΦ
        dδb = ck_aH*vb - 3*dΦ
        dvc = -vc - ck_aH*Ψ
        dvb = -vb - ck_aH*Ψ + τ′*R*(3*Θl[1]+vb)
        @assert lmax != 3 # TODO: how to calculate dΘl with lmax = 3?
        dΘl0   = -ck_aH*Θl[1] - dΦ
        dΘl1   =  ck_aH/3 * (Θ0-2*Θl[2]+Ψ) + τ′*(Θl[1]+vb/3)
        dΘl2   =  ck_aH/(2*2+1) * (2*Θl[2-1] - (2+1)*Θl[2+1]) + τ′*(Θl[2]-Π/10)
        dΘlmax = ck_aH*Θl[lmax-1] - (lmax+1)/(aH(co,x)*η(co,x))*Θl[lmax] + τ′*Θl[lmax] # 2nd term: their η is my c*η

        # re-pack variables into vector
        # TODO: pack and assign simultaneously?
        dy[i_δc]    = dδc
        dy[i_vc]    = dvc
        dy[i_δb]    = dδb
        dy[i_vb]    = dvb
        dy[i_Φ]     = dΦ
        dy[i_Θl(0)] = dΘl0
        dy[i_Θl(1)] = dΘl1
        dy[i_Θl(2)] = dΘl2
        for l in 3:lmax-1
            dy[i_Θl(l)] = ck_aH/(2*l+1) * (l*Θl[l-1] - (l+1)*Θl[l+1]) + τ′*(Θl[l]-0) # probably inlined? TODO: reduce allocation?
        end
        dy[i_Θl(lmax)] = dΘlmax
        return nothing # dy is in-place
    end

    # TODO: why do explicit methods take small steps when x ≈ -7.4 or so?
    # TODO: which quantity behaves weirdly here?
    # TODO: try to integrate with an explicit solver to a small time,
    # TODO: then plot all functions on a small interval around it to see if any of them behaves badly
    # TODO: specify stiff solver explicitly, or use an automatic one?
    y1 = initial_conditions_untight(co, x1, k, lmax)
    return _spline_integral(dy_dx!, x1, x2, y1; abstol=1e-9, reltol=1e-9, name="perturbations untight (k=$(k*Mpc)/Mpc)")
end

# TODO: join them
# TODO: if x_tight_latest = NaN, only integrate the full equations
function splined_perturbations_combined(co::ΛCDM, k::Real; lmax::Integer=30)
    spl1s = splined_perturbations_tight(co, k; lmax=lmax)
    spl2s = splined_perturbations_untight(co, k; lmax=lmax)
    spls = Vector{Spline1D}(undef, i_max(lmax))
    for i in 1:i_max(lmax)
        spls[i] = splinejoin(spl1s[i], spl2s[i])
    end
    return spls
end

# TODO: for tight and untight?
function splined_perturbations(co::ΛCDM; lmax::Integer=30)
    if length(co.perturbation_splines) == 0
        # then make it
        kmin, kmax = 0.00005 / Mpc, 0.3 / Mpc
        ks = range(kmin, kmax, 5) # TODO: what spacing?
        
        # take x values from most rapidly oscillating smallest-scale solution (k = kmax)
        perturbs_kmax = splined_perturbations_combined(co, kmax; lmax=lmax) # fill perturbations_untight_spline
        xs = unique(perturbs_kmax[1].t) # unique values only

        perturbs = Array{Float64, 3}(undef, length(perturbs_kmax), length(xs), length(ks)) # indexed as [i_quantity, i_x, i_k]
        for (i_k, k) in enumerate(ks)
            println("k = $(k*Mpc) / Mpc")
            p = splined_perturbations_combined(co, k; lmax=lmax) # fill perturbations_untight_spline
            for i_q in 1:length(perturbs_kmax)
                perturbs[i_q, :, i_k] .= p[i_q].(xs)
            end
        end

        for i_qty in 1:i_max(lmax)
            qty = perturbs[i_qty, :, :]
            push!(co.perturbation_splines, Spline2D(xs, ks, qty)) # spline (x, k)
        end
    end
    return co.perturbation_splines
end

δc(co::ΛCDM, x::Real, k::Real) = splined_perturbations(co)[i_δc](x, k)
δb(co::ΛCDM, x::Real, k::Real) = splined_perturbations(co)[i_δb](x, k)
vc(co::ΛCDM, x::Real, k::Real) = splined_perturbations(co)[i_vc](x, k)
vb(co::ΛCDM, x::Real, k::Real) = splined_perturbations(co)[i_vb](x, k)
 Φ(co::ΛCDM, x::Real, k::Real) = splined_perturbations(co)[i_Φ](x, k)
 Ψ(co::ΛCDM, x::Real, k::Real) = -Φ(co,x,k) - 12*co.H0^2/(c*k*a(x))^2 * (co.Ωγ0*Θl(co,x,k,2)) # TODO: neutrinos
Θl(co::ΛCDM, x::Real, k::Real, l::Integer) = splined_perturbations(co)[i_Θl(l)](x, k)
