include("Constants.jl")

# TODO: start/end with SI units, but solve equation in Planck units
# TODO: no, cannot do that! my other functions use SI units

# index map (1-based)
const i_δc = 1
const i_vc = 2
const i_δb = 3
const i_vb = 4
const i_Φ  = 5
const i_Θ0 = 6
i_Θl(l::Integer) = i_Θ0 + l
i_max(lmax::Integer) = i_Θl(lmax)

# TODO: include neutrinos
# TODO: include polarization
# TODO: check units
function initial_conditions(co::ΛCDM, x0::Real, k::Real, lmax::Integer)
    aH0 = aH(co, x0)
    dτ0 = dτ(co, x0)

    fν = 0 # TODO: include neutrinos: co.Ων0 / co.Ωr0
    Ψ  = -1 / (3/2 + 2*fν/5)
    Φ  = -(1 + 2*fν/5) * Ψ # TODO: "acts as normalization" ?
    δb = δc = -3/2 * Ψ
    vb = vc = -k*c/(2*aH0) * Ψ # TODO: units?
    function Θl(l::Integer)
        if l == 0
            return -1/2 * Ψ
        elseif l == 1
            return -c*k / (3*aH0) * Θl(0) # = k / (6*aH0) * Ψ
        elseif l == 2
            return -20*c*k / (45*aH0*dτ0) * Θl(1) # TODO: include polarization
        else
            return -l/(2*l+1) * c*k/(aH0*dτ0) * Θl(l-1)
        end
    end

    y = Vector{Float64}(undef, i_max(lmax))
    y[i_δc] = δc
    y[i_vc] = vc
    y[i_δb] = δb
    y[i_vb] = vb
    y[i_Φ]  = Φ
    for l in 0:lmax
        y[i_Θl(l)] = Θl(l)
    end

    return y
end

function time_tight_coupling(co::ΛCDM, k::Real)
    x1 = find_zero(x -> abs(dτ(co,x)) - 10,                (-20, +20))
    x2 = find_zero(x -> abs(dτ(co,x)) - 10 * c*k/aH(co,x), (-20, +20))
    x3 = -8.3 # TODO: time_recombination() or -8.3? at least a dynamic way of computing it?
    #x3 = time_recombination(co)
    #println("$x1 $x2 $x3")
    return min(x1, x2, x3)
end

# TODO: tight coupling is equivalent to lmax=2, then post-compute l>lmax
function perturbations_tight_coupling(co::ΛCDM, x::Real, k::Real; x0::Real=-20.0, lmax::Integer=1)
    if isnothing(co.perturb_tc_spline)
        # TODO: integrate dy/dx = f
        # input k is in SI-units (meters), but internal functions work with Planck units
        function dy_dx(x, y)
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
            Θ2 = -20/45 * ck_aH * Θ1/dτ(co,x)

            # 2) compute intermediate quantities
            Ψ = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θ2) # TODO: neutrinos

            # 3) compute derivatives
            dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0) # TODO: neutrinos
            dδc = ck_aH*vc - 3*dΦ
            dδb = ck_aH*vb - 3*dΦ
            dvc = -vc - ck_aH*Ψ
            dΘ0 = -ck_aH*Θ1 - dΦ
            #dΘ1 = ck_aH/3*(Θ0-2*Θ2+Ψ) + dτ(co,x) * (Θ1 + vb/3)
            q = (-((1-R)*dτ(co,x)+(1+R)*d2τ(co,x))*(3*Θ1+vb) - ck_aH*Ψ + (1-daH_aH)*ck_aH*(-Θ0+2*Θ2) - ck_aH*dΘ0) /
                ((1+R)*dτ(co,x) + daH_aH - 1)
            #dvb = -vb - ck/aH*Ψ + dτ(co,x) * R * (3*Θ1+vb)
            dvb = 1/(1+R) * (-vb - ck_aH*Ψ + R*(q+ck_aH*(-Θ0+2*Θ2)) - ck_aH*Ψ)
            dΘ1 = (q - dvb) / 3

            # re-pack variables into vector
            dy = Vector{Float64}(undef, i_max(lmax))
            dy[i_δc] = dδc
            dy[i_vc] = dvc
            dy[i_δb] = dδb
            dy[i_vb] = dvb
            dy[i_Φ]  = dΦ
            dy[i_Θl(0)] = dΘ0
            dy[i_Θl(1)] = dΘ1

            return dy
        end

        y0 = initial_conditions(co, x0, k, lmax)
        println("y0 = $y0")

        co.perturb_tc_spline = _spline_integral(dy_dx, x0, time_tight_coupling(co,k), y0)
    end

    return co.perturb_tc_spline(x)
end

# TODO: extend outside tight coupling
δc(co::ΛCDM, x::Real, k::Real)             = perturbations_tight_coupling(co, x, k)[i_δc]
δb(co::ΛCDM, x::Real, k::Real)             = perturbations_tight_coupling(co, x, k)[i_δb]
vb(co::ΛCDM, x::Real, k::Real)             = perturbations_tight_coupling(co, x, k)[i_vb]
vc(co::ΛCDM, x::Real, k::Real)             = perturbations_tight_coupling(co, x, k)[i_vc]
 Φ(co::ΛCDM, x::Real, k::Real)             = perturbations_tight_coupling(co, x, k)[i_Φ]
Θl(co::ΛCDM, x::Real, k::Real, l::Integer) = perturbations_tight_coupling(co, x, k)[i_Θl(l)]
