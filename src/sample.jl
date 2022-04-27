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
	for i in 1:L
		E += g.J[conf[i], a, i, j] - g.J[conf[i], b, i, j]
	end
	# Changing conf if needed
	if E <= 0. || exp(-E) > rand(rng)
		conf[j] = b
	end
	return nothing
end

