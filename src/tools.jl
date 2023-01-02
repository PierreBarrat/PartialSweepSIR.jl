"""
	cross_immunity(N::Int, variant_list)

Cross-immunity matrix. Each element in `variant_list` is interpreted as a pair of the form
`(b,f)` giving the cross-immunity of this variant to all others.
The value `variant_list[1]` for the first variant is not used: there are no previous
variants, making them useless.

### Example

```
cross_immunity(2, [(b=0.8, f=0.5)])
cross_immunity(2, [(0., 0.), (b=0.8, f=0.5)]) # first value of variant_list ignored
cross_immunity(3, [(0.6, 0.5), (0.3, 0.2)])
"""
function cross_immunity(N::Int, variant_list)
	K = diagm(ones(N))
	# custom function for indexing into variant_list
	idx(b) = length(variant_list) == N ? b : b-1
	for b in 1:N, a in 1:(b-1)
		# b > a --> a ancient & b recent
		β, ϕ = variant_list[idx(b)]
		K[a,b] = β # immunity from b infection to variant a
		K[b,a] = ϕ # immunity from a infection to variant b
	end
	return K
end

"""
	cross_immunity(N::Int, backward::Real, forward::Real)

Cross-immunity matrix with the same backward and forward effect for all variants.
"""
function cross_immunity(N::Int, backward::Real, forward::Real)
	K = diagm(ones(N))
	for a in 1:N, b in (a+1):N
		K[a,b] = backward # cross-immunity from next virus to previous virus (a < b)
		K[b,a] = forward
	end
	return K
end

"""
	cross_immunity(N)

Diagonal cross-immunity matrix (*i.e.* no interaction between variants).
"""
cross_immunity(N::Int) = diagm(ones(N))

"""
	grow_cross_immunity(K, b, f)

Return a matrix of dimension (N+1)x(N+1), with `N=size(K,1)`

K | b
-----
f | 1

"""
function grow_cross_immunity(K, b, f)
	N = size(K,1)
	new_K = zeros(N+1, N+1)
	new_K[1:N, 1:N] .= K
	new_K[N+1, 1:N] .= b
	new_K[1:N, N+1] .= f
	new_K[N+1, N+1] = 1
	return new_K
end



"""
	frequency(X::SIRState, region::Int, virus::Int)

Fraction of infections caused by `virus` in `region`.
"""
function frequency(X::SIRState, region, virus)
	return _frequency(X, region, virus)
end
"""
	frequency(X::SIRState, virus::Int)

Frequency of `virus` accross regions.
"""
function frequency(X::SIRState, virus::Int)
	return mean(frequency(X, i, virus) for i in 1:length(regions(X)))
end

_frequency(X::SIRState, i::Int, a) = X.regions[i].I[a] / sum(X.regions[i].I)
function _frequency(X::SIRState, i::Colon, a)
	return [_frequency(X, ii, a) for ii in 1:length(regions(X))]
end

"""
	frequency(sol::Solution, tvals, i, a)

Frequency of infections by virus `a` in region `i`.
"""
function frequency(sol::Solution, tvals, i, a)
	I = sol[tvals, i, :, :I]
	return map(x -> x[a], I) ./ map(sum, I)
end
"""
	frequency(sol::Solution, tvals, a)

Frequency of infections by virus `a` accross regions.
"""
function frequency(sol::Solution, tvals, a)
	mean(frequency(sol, tvals, i, a) for i in 1:sol.parameters.M)
end





"""
	set_infected(r::Region, a::Int, val)

Set `r.viruses[a].I` to `val`. Equilibrate host pop by taking from `S`.
Useful to introduce a new virus.
"""
function _set_infected!(r::Region, a::Int, val)
	old_v = r.viruses[a]
	new_v = Virus(;
		I = val,
		C = old_v.C,
		R = old_v.R,
		S = 1 - old_v.C - old_v.R - val
	)

	r.viruses[a] = new_v

	return nothing
end
"""
	set_infected(X::SIRState, a::Int, val)

Return a copy of `X` with `viruses[a].I` set to `val` in all regions.
"""
function set_infected(X::SIRState, a::Int, val)
	X_copy = deepcopy(X)
	for r in regions(X_copy)
		_set_infected!(r, a, val)
	end
	return X_copy
end
