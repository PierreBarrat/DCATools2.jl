"""
	rand_pottsgraph(L, q; rJ = randn, rh = randn)

Return random `PottsGraph` of dimensions `(L, q)`. Kwargs `rJ` and `rh` define the random function used for elements of `J` and `h`

## Examples

A random graph of dimensions `L=3` and `q=2` with normally distributed parameters, with mean/standard deviation `(0.5, 0.1)` for fields and `(0., 0.01)` for couplings.
```
rand_pottsgraph(
	3, 2; 
	rJ = ()->0.01*randn(), rh = ()->0.1*randn() + 0.5
)
```
"""
function rand_pottsgraph(L, q; rJ = randn, rh = randn)
	g = PottsGraph(L, q)
	for i in 1:L, j in (i+1):L, a in 1:q, b in 1:q
		g.J[a, b, i, j] = rJ()
		g.J[b, a, j, i] = g.J[a, b, i, j]
	end
	for i in 1:L
		g.J[:, :, i, i] .= 0
	end
	for i in 1:L, a in 1:q
		g.h[a, i] = rh()
	end

	return g
end

"""
	matrix_representation(g::PottsGraph)

Return the matrix/vector representation of the parameters of `g`: couplings are represented
by a `(L*q) x (L*q)` matrix and fields by a `L*q` vector.
"""
function matrix_representation(g::PottsGraph)
	q, L = size(g)
	J = zeros(Float64, L*q, L*q)
	h = zeros(Float64, L*q)

	for i in 1:L, j in 1:L, a in 1:q, b in 1:q
		J[(i-1)*q + a, (j-1)*q + b] = g.J[a,b,i,j]
	end
	for i in 1:L, a in 1:q
		h[(i-1)*q + a] = g.h[a, i]
	end

	return J, h
end