module Algorithms

export MetropolisHastings

using Distributions

# TODO: burnin?
# TODO: constant step scaling "locks" their mutual *ratios*; scale only one random step?

# TODO: wtf have i done now?
function MetropolisHastings(logL::Function, bounds::Tuple{Vector{Float64},Vector{Float64}}, samples::Integer; steps=nothing, verbose=true)
    params_lo, params_hi = bounds

    @assert !isnan(logL(params_lo)) && !isnan(logL(params_hi))
    @assert all(params_hi .>= params_lo)
    @assert length(params_lo) == length(params_hi)
    @assert steps == nothing || length(steps) == length(params_lo)

    # pre-allocate arrays to record accepted parameters and likelihoods
    nparams = length(params_lo)
    params = Array{Float64, 2}(undef, samples, nparams)
    logLs = Vector{Float64}(undef, samples)

    # unless step sizes are explicitly specified, automatically set them to 5% of bounds
    if steps == nothing
        steps = 0.05 .* (params_hi .- params_lo)
    end

    ndistr = Normal()
    udistr = Uniform(0, 1)

    # set first params to a random guess within the parameter bounds
    curr_params = params_lo .+ rand.(udistr) .* (params_hi .- params_lo)
    curr_logL = logL(curr_params)
    if verbose
        println("Step sizes: ", steps)
        print("Initial random guess: ")
        println("params = ", curr_params, ", log(L) ∝ ", round(curr_logL, digits=2))
    end

    step = 0
    sample = 0
    while sample < samples
        step += 1
        #new_params .= curr_params .+ rand.(ndistr) .* steps
        new_params = curr_params .+ rand.(ndistr) .* steps # TODO: optimize
        #println(new_params)
        inbounds = all(i -> params_lo[i] <= new_params[i] <= params_hi[i], 1:nparams)
        ##inbounds = all(params_lo .<= new_params .<= params_hi)
        new_logL = inbounds ? logL(new_params) : -Inf

        # accept new state unconditionally if new likelihood is greater; otherwise with probability new_L/curr_L
        if new_logL - curr_logL > log(rand(udistr)) # equivalent to new_L > rand(udistr) * curr_L (TODO: optimize: avoid rand() every loop)
            sample += 1
            curr_params = new_params
            curr_logL = new_logL

            params[sample,:] = curr_params # record params
            logLs[sample] = curr_logL # record likelihood

            if verbose
                accrate = sample / step
                print("\33[2K\r") # overwrite current line
                print("Sample ", sample, "/", samples, " (accept ", round(accrate*100, digits=1), "%): ")
                print("params = ", curr_params, ", log(L) ∝ ", round(curr_logL, digits=2))
                flush(stdout) # flush to display when line does not end with newline
            end
        end

        # After running for a while, check if we are accepting too many or few states
        # If so, tune step sizes and re-run: https://colcarroll.github.io/hmc_tuning_talk/
        # This *re-runs the whole algorithm*, so it *does not ruin the Markov chain nature*, but is merely equivalent to fiddling with step sizes
        if step >= samples # have run for some time
            accrate = sample / step
            if abs(accrate - 0.25) > 0.10 # aim for acceptance rate 25%: https://doi.org/10.1214/aoap/1034625254
                scale = 2^((accrate-0.25)/(0.50-0.25)) # 0.5x/2x when accrate is 0%/50% (TODO: sensible?)
                steps .*= scale
                if verbose
                    println("\nToo far from 25% accept rate; re-running with step sizes scaled by ", scale)
                end
                return MetropolisHastings(logL, bounds, samples; steps, verbose)
            end
        end
    end

    if verbose
        println() # end the one line that we are constantly overwriting
    end

    return params, logLs
end

end
