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
println("($(length(pspec.k)) k-values)")

for Δx in [0.01, 0.05]
    print("CMB power spectrum (Δx = $(Δx)):")
    pspec = @btime CMBPowerSpectrum(rec, :TT; xs=range(-10,0,step=$Δx))
    println("($(length(pspec.l))) l-values")
end

end
