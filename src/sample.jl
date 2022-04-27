"""
	mcmc_step!(conf, g[, rng=GLOBAL_RNG])

Perform one MCMC flips to `conf`.
"""
function mcmc_step!(conf, g, rng=Random.GLOBAL_RNG)
	q, L = size(g)
	E = 0.
	# Flip position
	j = rand(rng, 1:L)
	a = conf[j] # initial state
	b = mod(a - 1 + rand(rng, 1:(q-1)), q) + 1 # new state
	# Computing energy change - remember E = - (J + h)
	E += g.h[a, j] - g.h[b, j]
	@inbounds for i in 1:L
		# E += g.J[conf[i], a, i, j]  - g.J[conf[i], b, i, j]
		E += g.J[a, conf[i], j, i] - g.J[b, conf[i], j, i]
	end
	# Changing conf if needed
	if E <= 0. || exp(-E) > rand(rng)
		conf[j] = b
	end
	return nothing
end

"""
	mcmc_sweep!(conf, g[, rng = Random.GLOBAL_RNG])

Perform `L` MCMC flips to `conf` where `L` is the length of `conf`.
"""
function mcmc_sweep!(conf, g, rng = Random.GLOBAL_RNG)
	q, L = size(g)
	@assert length(conf) == L "Configuration and graph have incompatible sizes\
		(respectively $(length(conf)) and $L"
	for rep in 1:L
		mcmc_step!(conf, g, rng)
	end
	return nothing
end

"""
	sample(
		g, M; 
		init = rand(1:size(g)[1], size(g)[2]), 
		Twait = 1, 
		burnin = 5*Twait, 
		rng = Random.GLOBAL_RNG
	)

Sample `M` configurations from `g` using MCMC. 

## Kwargs
- `init`: initial configuration
- `Twait`: number of sweeps between two sample, *i.e.* decorrelation time of samples
- `burnin`: number of sweeps between initial configuration and first sample
"""
function sample(
	g, M; 
	init = rand(1:size(g)[1], size(g)[2]), 
	Twait = 1, 
	burnin = 5*Twait, 
	rng = Random.GLOBAL_RNG
)
	q, L = size(g)
	@assert length(init) == L "Init. configuration and graph have incompatible sizes\
		(respectively $(length(init)) and $L"
	conf = copy(init)
	sample = zeros(Int64, M, L)

	for t in 1:burnin
		mcmc_sweep!(conf, g, rng)
	end
	sample[1, :] .= conf
	for it in 2:M
		for t in 1:Twait
			mcmc_sweep!(conf, g, rng)
		end
		sample[it, :] .= conf
	end

	return sample
end
















