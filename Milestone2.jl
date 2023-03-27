module Milestone2

include("Cosmology.jl")
using .Cosmology

co = ΛCDM()

println(Tγ(co, 0))
println(Xe_saha(co, -10))
println(Xe_saha(co, 0))

end
