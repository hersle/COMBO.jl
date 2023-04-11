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
            return -l/(2*l+1) * k/(aH0*dτ0) * Θl(l-1)
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

function perturbations(co::ΛCDM, k::Real, x0::Real, lmax::Integer)
    # TODO: integrate dy/dx = f
    # input k is in SI-units (meters), but internal functions work with Planck units
    y0 = initial_conditions(co, k, x0, lmax)
    println("y0 = $y0")
end
