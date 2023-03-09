module Milestone2

include("Cosmology.jl")
using .Cosmology

co = ΛCDM()

println("hi")
println(Tγ(co, 0))

end
