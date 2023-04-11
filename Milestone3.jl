module Milestone3

include("Cosmology.jl")
include("Constants.jl")
include("Plotting.jl")

using .Cosmology
using .Constants
using LaTeXStrings

co = ΛCDM()
perturbations(co, 0.1 / Mpc, 1e-18, 10)

end
