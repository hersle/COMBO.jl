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
function initial_conditions(co::ΛCDM, k::Real, x0::Real, lmax::Integer)
    aH0 = aH(co, x0)
    dτ0 = dτ(co, x0)

    fν = 0 # TODO: include neutrinos: co.Ων0 / co.Ωr0
    Ψ  = -c^2 / (3/2 + 2*fν/5)
    Φ  = -(1 + 2*fν/5) * Ψ # TODO: "acts as normalization" ?
    δb = δc = -3/2 * Ψ / c^2
    vb = vc = -k/(2*aH0) * Ψ # TODO: units?

    Θ0 = -1/2 * Ψ
    Θ1 = -k / (3*aH0) * Θ0 # = k / (6*aH0) * Ψ
    Θ2 = -20*k / (45*aH0*dτ0) * Θ1 # TODO: include polarization

    y = Vector{Float64}(undef, i_max(lmax))
    y[i_δc] = δc
    y[i_vc] = vc
    y[i_δb] = δb
    y[i_vb] = vb
    y[i_Φ]  = Φ
    y[i_Θl(0)] = -1/2 * Ψ
    y[i_Θl(1)] = -k / (3*aH0) * Θ0 # = k / (6*aH0) * Ψ
    y[i_Θl(2)] = -20*k / (45*aH0*dτ0) * Θ1 # TODO: include polarization
    for l in 3:lmax
        y[i_Θl(l)] = -l/(2*l+1) * k/(aH0*dτ0) * y[i_Θl(l-1)]
    end

    return y
end

function perturbations(co::ΛCDM, k::Real, x0::Real, lmax::Integer)
    # TODO: integrate dy/dx = f
    # input k is in SI-units (meters), but internal functions work with Planck units
    y0 = initial_conditions(co, k, x0, lmax)
    println("y0 = $y0")
end
