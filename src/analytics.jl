"""
	equilibrium(X::SIRState)

Return an `SIRState` object giving the equilibrium reached by the system if initiated
from `X`.
"""
function equilibrium(X::SIRState)
	return if length(regions(X)) == 1
		equilibrium_1(X)
	else
		equilibrium_M(X)
	end
end

function equilibrium_1(X::SIRState)
	@assert X.parameters.M == 1 "Only for one region"
	region = _equilibrium(X.regions[1], X.parameters)
	return SIRState(; regions = [region], parameters = deepcopy(X.parameters))
end

function _equilibrium(r::Region, p::Parameters)
	S, I, C, R = _equilibrium(r.K, p)
	viruses = Vector{Virus}(undef, p.N)
	for (a, (s, i, c, r)) in enumerate(zip(S, I, C, R))
		viruses[a] = Virus(; S=s, I=i, C=c, R=r)
	end
	return Region(; viruses, K=copy(r.K))
end
function _equilibrium(K::Matrix, p::Parameters)
	@unpack N, α, γ, δ = p

	S = δ/α * ones(N)
	I = try
		γ/δ * (1-δ/α) * (K \ ones(N))
	catch err
		@error "Singular `K`: rank = $(rank(K)) < $N. Cannot compute equilibrium." K
		error(err)
	end
	C = γ/δ * (1-δ/α) .- I
	R = 1 .- S .- I .- C

	return S, I, C, R
end





# """
# !! THIS IS FALSE FOR S and R (but should be good for I)
# """
function equilibrium_M(X::SIRState)
	@unpack M, N, α, γ, δ, C = X.parameters

	eq_regions = []
	for (i, region) in enumerate(X.regions)
		Z = sum(C[i,:])
		K_av = sum(C[i,j] * r.K for (j,r) in enumerate(X.regions)) / Z

		S = δ/α * ones(N)
		I = try
			γ/δ * (1-δ/α) * (K_av \ ones(N))
		catch err
			@error "Singular `K_av`: rank = $(rank(K_av)) < $N. Cannot compute equilibrium." K_av
			error(err)
		end
		R = 1 .- S .- I

		viruses = Vector{Virus}(undef, N)
		for (a, (s, i, c, r)) in enumerate(zip(S, I, C, R))
			viruses[a] = Virus(; S=s, I=i, C=0, R=r)
		end

		push!(eq_regions, Region(; viruses, K=copy(region.K)))
	end

	return SIRState(; regions = eq_regions, parameters = deepcopy(X.parameters))
end


