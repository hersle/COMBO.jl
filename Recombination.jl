include("Constants.jl")

ρcrit(co::ΛCDM, x::Real) = 3 * H(co, x)^2 / (8 * π * G)
ρb(co::ΛCDM, x::Real) = Ωb(co, x) * ρcrit(co, x)
nb(co::ΛCDM, x::Real) = ρb(co, x) / mH
Tγ(co::ΛCDM, x::Real) = co.Tγ0 / a(x)

function Xe_saha(co::ΛCDM, x::Real)
    # TODO: add He
    # TODO: are there any overflow issues? plot looks ok.
    T = Tγ(co, x)
    a = 1
    b = (2 * π * me * kB * T / h^2)^(3/2) * exp(-EHion/(kB*T)) / nb(co,x)
    c = -b
    return quadroots(a, b, c)[2]
end

function Xe_peebles(co::ΛCDM, x::Real)
    # TODO:
end

function Xe(co::ΛCDM, x::Real)
    # TODO: adaptively branch into saha/peebles
end
