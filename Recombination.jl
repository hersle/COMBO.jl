include("Constants.jl")

ρcrit(co::ΛCDM, x::Real) = 3 * H(co, x)^2 / (8 * π * G)
ρb(co::ΛCDM, x::Real) = Ωb(co, x) * ρcrit(co, x)
nb(co::ΛCDM, x::Real) = ρb(co, x) / mH
Tγ(co::ΛCDM, x::Real) = co.Tγ0 / a(x)

function xe_saha(co::ΛCDM, x::Real)
    T = Tγ(co, x)
    a = 1
    b = (2 * π * me * kB * T / h^2)^(3/2) * exp(-EHion/(kB*T)) / nb(co,x)
    c = -b

    d = b^2 - 4*a*c # >= b^2, since a > 0 and c < 0!
    @assert b >= 0

    return (-b + √(d)) / (2*a) # - solution is always negative
end
