module MilestoneBenchmark

include("Cosmology.jl")
include("Constants.jl")
using .Cosmology
using .Constants
using BenchmarkTools


print("Parameters:")
par = @btime Parameters()

print("Background:")
bg = @btime Background(par)

print("Recombination:")
rec = @btime Recombination(bg)

for k in [0.01, 1, 100] / Mpc
    print("Perturbations (k = $(k*Mpc)/Mpc):")
    perturb = @btime PerturbationMode(rec, $k; verbose=false)
end

print("Matter power spectrum:")
pspec = @btime MatterPowerSpectrum(rec)

print("CMB power spectrum:")
pspec = @btime CMBPowerSpectrum(rec, :TT)

end
