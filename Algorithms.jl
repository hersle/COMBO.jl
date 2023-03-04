module Algorithms

export MetropolisHastings

using Distributions

# TODO: constant step scaling "locks" their mutual *ratios*; scale only one random step? use covariance matrix somehow?
function MetropolisHastings(logL::Function, bounds::Vector{Tuple{Float64,Float64}}, samples::Integer, chains::Integer; steps=nothing, burnin::Integer=0, verbose=true)
    params_lo = [bound[1] for bound in bounds]
    params_hi = [bound[2] for bound in bounds]

    @assert !isnan(logL(params_lo)) && !isnan(logL(params_hi))
    @assert all(params_hi .>= params_lo)
    @assert length(params_lo) == length(params_hi)
    @assert steps == nothing || length(steps) == length(params_lo)
    @assert chains > 0

    # pre-allocate arrays to record accepted parameters and likelihoods
    nparams = length(params_lo)
    params = Array{Float64, 2}(undef, chains * samples, nparams)
    logLs = Vector{Float64}(undef, chains * samples)

    if chains > 1 # run individual chains and concatenate their outputs
        for chain in 1:chains
            if verbose
                println("Metropolis-Hastings chain $chain/$chains")
            end
            params1, logLs1 = MetropolisHastings(logL, bounds, samples, 1; steps, burnin, verbose) # run one chain
            params[1+(chain-1)*samples:chain*samples, :] = params1
            logLs[1+(chain-1)*samples:chain*samples] = logLs1
        end
        return params, logLs
    end # otherwise, continue and just run one chain

    # unless step sizes are explicitly specified, automatically set them to 1% of bounds
    if steps == nothing
        steps = 0.01 .* (params_hi .- params_lo)
    end

    ndistr = Normal(0, 1)
    udistr = Uniform(0, 1)

    # set first params to a random guess within the parameter bounds
    curr_params = params_lo .+ rand.(udistr) .* (params_hi .- params_lo)
    curr_logL = logL(curr_params)

    if verbose
        println("Step sizes: ", steps)
        print("Initial guess: ")
        println("params = ", curr_params, ", log(L) ∝ ", round(curr_logL, digits=2))
    end

    step = 0
    sample = 0
    while sample < samples + burnin
        step += 1
        new_params = curr_params .+ rand.(ndistr) .* steps
        inbounds = all(i -> params_lo[i] <= new_params[i] <= params_hi[i], 1:nparams)
        new_logL = inbounds ? logL(new_params) : -Inf
        accept = new_logL - curr_logL > log(rand(udistr)) # TODO: why the fuck does this line allocate memory?
        if accept # equivalent to new_L/curr_L > rand(udistr) (TODO: avoid individual calls to rand() every loop?)
            sample += 1
            curr_params .= new_params
            curr_logL = new_logL

            # record sample parameters and likelihood, if past the burn-in stage
            if sample > burnin
                params[sample-burnin,:] .= curr_params
                logLs[sample-burnin] = curr_logL
            end
        end

        accrate = sample / step

        if verbose && accept
            print("\33[2K\r") # overwrite current line
            print("Sample ", sample-burnin, "/", samples) # sample-burnin < 0 indicates burn-in stage
            print(" (accept ", round(accrate*100, digits=1), "%): ")
            print("params = ", curr_params, ", log(L) ∝ ", round(curr_logL, digits=2))
            flush(stdout) # flush to display when line does not end with newline
        end

        # After running for a while, check if we are accepting too many or few states
        # If so, tune step sizes and re-run: https://colcarroll.github.io/hmc_tuning_talk/
        # This *re-runs the whole algorithm*, so it *does not ruin the Markov chain nature*
        # Aim for acceptance rate 25%: https://doi.org/10.1214/aoap/1034625254
        if step >= samples && abs(accrate - 0.25) > 0.10 # have run for some time and accept rate is far from target
            scale = 2^((accrate-0.25)/(0.50-0.25)) # 0.5x/2x when accrate is 0%/50% (TODO: sensible?)
            steps .*= scale
            if verbose
                println("\nToo far from 25% accept rate; re-running with step sizes scaled by ", scale)
            end
            return MetropolisHastings(logL, bounds, samples, 1; steps, burnin, verbose)
        end
    end

    if verbose
        println() # end that line that we are constantly overwriting
    end

    return params, logLs
end

end
