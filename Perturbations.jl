include("Constants.jl")

# index map (1-based)
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
    x3 = -8.3 # TODO: time_recombination() or -8.3? at least a dynamic way of computing it?
    return min(x1, x2, x3)
end

# tight coupling is equivalent to lmax=2, then post-compute l>lmax
# TODO: add tight_coupling flag, recurse to find full solution?
# TODO: or simply have an if-else in dy_dx to switch "on the fly"
# TODO: store splines for each requested k-value
function perturbations_tight(co::ΛCDM, x::Real, k::Real; x1::Real=-20.0)
    x2 = time_tight_coupling(co, k)
    @assert x1 <= x <= x2 "x = $x, x1 = $x1, x2 = $x2"

    if isnothing(co.perturbations_tight_spline)
        # integrate dy/dx = f
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

            # 2) compute derivatives
            Ψ   = -Φ - 12 * (co.H0 / (c*k*a(x)))^2 * (co.Ωγ0*Θ2) # TODO: neutrinos
            dΦ  = Ψ - ck_aH^2/3*Φ + (co.H0/aH(co,x))^2/2 * (co.Ωc0/a(x)*δc + co.Ωb0/a(x)*δb + 4*co.Ωγ0/a(x)^2*Θ0) # TODO: neutrinos
            dδc = ck_aH*vc - 3*dΦ
            dδb = ck_aH*vb - 3*dΦ
            dvc = -vc - ck_aH*Ψ
            dΘ0 = -ck_aH*Θ1 - dΦ
            q   = (-((1-R)*dτ(co,x)+(1+R)*d2τ(co,x))*(3*Θ1+vb) - ck_aH*Ψ + (1-daH_aH)*ck_aH*(-Θ0+2*Θ2) - ck_aH*dΘ0) /
                  ((1+R)*dτ(co,x) + daH_aH - 1)
            dvb = 1/(1+R) * (-vb - ck_aH*Ψ + R*(q+ck_aH*(-Θ0+2*Θ2)) - ck_aH*Ψ)
            dΘ1 = (q - dvb) / 3

            # re-pack variables into vector
            lmax = 1 # tight coupling: only need to integrate Θ(l<=2) (see above comment)
            dy = Vector{Float64}(undef, i_max(lmax))
            dy[i_δc]    = dδc
            dy[i_vc]    = dvc
            dy[i_δb]    = dδb
            dy[i_vb]    = dvb
            dy[i_Φ]     = dΦ
            dy[i_Θl(0)] = dΘ0
            dy[i_Θl(1)] = dΘ1
            return dy
        end

        y1 = initial_conditions_tight(co, x1, k)
        co.perturbations_tight_spline = _spline_integral(dy_dx, x1, x2, y1)
    end

    return co.perturbations_tight_spline(x)
end

# TODO: extend outside tight coupling
δc(co::ΛCDM, x::Real, k::Real)             = perturbations_tight(co, x, k)[i_δc]
δb(co::ΛCDM, x::Real, k::Real)             = perturbations_tight(co, x, k)[i_δb]
vb(co::ΛCDM, x::Real, k::Real)             = perturbations_tight(co, x, k)[i_vb]
vc(co::ΛCDM, x::Real, k::Real)             = perturbations_tight(co, x, k)[i_vc]
 Φ(co::ΛCDM, x::Real, k::Real)             = perturbations_tight(co, x, k)[i_Φ]
Θl(co::ΛCDM, x::Real, k::Real, l::Integer) = perturbations_tight(co, x, k)[i_Θl(l)]
