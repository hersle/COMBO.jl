module Algorithms

export MetropolisHastings

using Distributions

# TODO: adaptive step size?
function MetropolisHastings(logL::Function, p_lo::Vector{Float64}, p_hi::Vector{Float64}; samples::Int=1000, steps = nothing, callback = (step, sample, samples, params) -> ())
    @assert !isnan(logL(p_lo)) && !isnan(logL(p_hi))
    @assert all(p_hi .>= p_lo)
    @assert length(p_lo) == length(p_hi)
    @assert steps == nothing || length(steps) == length(p_lo)

    nparams = length(p_lo)
    params = Array{Float64, 2}(undef, samples, nparams)
    logLs = Vector{Float64}(undef, samples)

    # TODO: seed here?
    ndist = Normal()
    udist = Uniform(0, 1)

    # random initial guess within bounds
    curr_params = p_lo .+ rand.(udist) .* (p_hi .- p_lo)
    curr_logL = logL(curr_params)
    println("Initial guess: $curr_params, logL = $curr_logL")

    # adaptive step sizes?
    steps = steps == nothing ? 0.05 .* (p_hi .- p_lo) : steps # TODO: how to choose scale?

    step = 0
    sample = 0
    while sample < samples
        step += 1
        new_params = curr_params .+ rand.(ndist) .* steps # TODO: handle initial guess nicely?
        new_logL = all(p_lo .<= new_params .<= p_hi) ? logL(new_params) : -Inf

        # new_L > curr_L: accept new state unconditionally (since rand(udist) < 1)
        # new_L < curr_L: accept new state with probability new_L / curr_L
        if new_logL - curr_logL > log(rand(udist)) # if new_L > rand(udist) * curr_L
            sample += 1
            curr_params = new_params
            curr_logL = new_logL

            callback(step, sample, samples, curr_params) # execute any user-provided callback
            params[sample,:] = curr_params # record params
            logLs[sample] = curr_logL # record likelihood
        end

        # After we have run for a while, check if we are accepting too many or few states
        # and re-run with different step sizes, if so: https://colcarroll.github.io/hmc_tuning_talk/
        # This *re-runs the whole algorithm*, so it *does not ruin the Markov chain nature*
        # It is merely equivalent to manually fiddling with the step sizes
        # TODO: the constant scaling "locks" the *ratio* between the parameters
        # TODO: select one random parameter and scale it only?
        if step == 3 * samples # have run for some time (collected 10% of samples)
            accrate = sample / step
            if abs(accrate - 0.25) > 0.10 # aim for acceptance rate 25%: https://doi.org/10.1214/aoap/1034625254
                scale = 2^((accrate-0.25)/(0.50-0.25))
                steps = scale .* steps # 0.5x/2x when accrate is 0%/50% (TODO: sensible?)
                println("Multiplying stepsize by $scale")
                return MetropolisHastings(logL, p_lo, p_hi; samples, steps, callback)
            end
        end
    end

    return params, logLs
end

end
